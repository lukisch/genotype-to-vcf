#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import gzip
import json
import time
import random
import string
import shutil
import socket
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import requests
from pyfaidx import Faidx
import psutil

# PyQt6 Imports
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QComboBox, 
                             QTextEdit, QProgressBar, QFileDialog, QMessageBox, 
                             QFrame, QStyleFactory)
from PyQt6.QtCore import QThread, pyqtSignal, Qt, QObject
from PyQt6.QtGui import QColor, QPalette, QFont, QIcon

# -----------------------------
# Konfiguration
# -----------------------------
APP_TITLE = "Genotype → VCF Pro Converter"
APP_VERSION = "1.0.0"
FASTA_URLS = {
    "GRCh37": "http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
    "GRCh38": "http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
}
FASTA_PATHS = {
    "GRCh37": "Homo_sapiens.GRCh37.dna.primary_assembly.fa",
    "GRCh38": "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
}
DBSNP_URL = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
TARGET_CPU_FRACTION = 0.70
MIN_WORKERS, MAX_WORKERS = 4, 200
TIMEOUT_SEC, RETRY_MAX = 15, 4
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.path.join(BASE_DIR, "cache.json")

_cache_lock = threading.Lock()

BUILD_HINTS_37 = [r"\bgrch\s*37\b", r"\bgrch37\b", r"\bbuild\s*37\b", r"\bhg19\b", r"\bb37\b", r"\bncb[iy]\s*build\s*37\b"]
BUILD_HINTS_38 = [r"\bgrch\s*38\b", r"\bgrch38\b", r"\bbuild\s*38\b", r"\bhg38\b", r"\bb38\b", r"\bncb[iy]\s*build\s*38\b"]

# -----------------------------
# Hilfsfunktionen & Logik
# -----------------------------
def is_rs_id(id_str: str) -> bool:
    return id_str.lower().startswith("rs") and id_str[2:].isdigit()

def atomic_write_json(path, obj):
    tmp = f"{path}.tmp.{''.join(random.choices(string.ascii_lowercase, k=6))}"
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(obj, f)
        for attempt in range(3):
            try:
                os.replace(tmp, path)
                break
            except PermissionError:
                if attempt == 2: raise
                time.sleep(0.5)
    except Exception as e:
        print(f"[ERROR] Fehler beim Schreiben von {path}: {e}")
        if os.path.exists(tmp): os.remove(tmp)

def cpu_count():
    try: return os.cpu_count() or 4
    except: return 4

def get_cpu_usage_percent():
    try: return psutil.cpu_percent(interval=0.1)
    except: return 30.0

def backoff_time(attempt):
    return min(30.0, (2 ** attempt) + random.random())

def load_cache(path=CACHE_FILE):
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except: pass
    empty = {}
    atomic_write_json(path, empty)
    return empty

def save_cache(cache, path=CACHE_FILE):
    atomic_write_json(path, cache)

def cache_upsert(cache, rsid, build, chrom, pos, ref):
    with _cache_lock:
        entry = cache.get(rsid) or {"assemblies": {}, "last_update": 0}
        entry["assemblies"][build] = {"chrom": chrom, "pos": int(pos), "ref": str(ref).upper()}
        entry["last_update"] = int(time.time())
        cache[rsid] = entry

def lookup_rsid_from_cache(chrom, pos, build, cache):
    chrom = str(chrom)
    pos = int(pos)
    for rid, entry in cache.items():
        hit = entry.get("assemblies", {}).get(build)
        if hit and str(hit["chrom"]) == chrom and int(hit["pos"]) == pos:
            return rid
    return None

# -----------------------------
# Genetik
# -----------------------------
def detect_sex_from_variants(variants):
    valid_y = 0
    threshold = 5
    for _, chrom, _, genotype in variants:
        chrom = chrom.upper().replace("CHR", "")
        if chrom == "Y":
            clean_gt = genotype.replace("_", "-").strip()
            if clean_gt not in ("--", "-", "00"):
                valid_y += 1
                if valid_y > threshold: return "male"
    return "female"

def in_par(chrom: str, pos: int, build: str) -> bool:
    chrom = chrom.upper().replace("CHR", "")
    if build == "GRCh37":
        par_regions = {"X": [(60001, 2699520), (154931044, 155260560)], "Y": [(10001, 2649520), (59034051, 59363566)]}
    elif build == "GRCh38":
        par_regions = {"X": [(10001, 2781479), (155701383, 156030895)], "Y": [(10001, 2781479), (56887903, 57217415)]}
    else:
        return False
    return any(start <= pos <= end for (start, end) in par_regions.get(chrom, []))

def ploidy_for_site(chrom: str, pos: int, build: str, sex: str) -> int:
    c = chrom.upper().replace("CHR", "")
    s = (sex or "unknown").lower()
    if c in {str(i) for i in range(1, 23)}: return 2
    if c == "MT" or c == "M": return 1
    if c == "X":
        if s == "female": return 2
        if s == "male": return 2 if in_par(c, pos, build) else 1
        return 2
    if c == "Y":
        if s == "female": return 0
        if s == "male": return 2 if in_par(c, pos, build) else 1
        return 0
    return 2

def gt_from_bases_snp(alleles, ref_base: str, alt_list, ploid: int):
    allele_map = {ref_base: "0"}
    for i, a in enumerate(alt_list, start=1):
        allele_map[a] = str(i)
    if ploid == 2:
        if len(alleles) != 2: return None
        gt_parts = [allele_map.get(a, ".") for a in alleles]
        gt = "/".join(gt_parts)
        return None if "." in gt else gt
    elif ploid == 1:
        if len(alleles) < 1: return None
        return allele_map.get(alleles[0], None)
    return None

# -----------------------------
# FASTA
# -----------------------------
def build_fasta_index(fasta_path):
    if not os.path.exists(fasta_path + ".fai"):
        Faidx(fasta_path)
    return fasta_path + ".fai"

def load_fai_index(fai_path):
    idx = {}
    if os.path.exists(fai_path):
        with open(fai_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5: continue
                idx[parts[0]] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
    return idx

def fetch_base_from_fasta(fasta_path, fai_index, chrom, pos):
    if chrom not in fai_index:
        if "chr" + chrom in fai_index: chrom = "chr" + chrom
        else: return "N"
    length, offset, line_bases, line_width = fai_index[chrom]
    if pos < 1 or pos > length: return "N"
    pos0 = pos - 1
    newlines = pos0 // line_bases
    byte_pos = offset + pos0 + (newlines * (line_width - line_bases))
    try:
        with open(fasta_path, "rb") as f:
            f.seek(byte_pos)
            return f.read(1).decode("ascii").upper()
    except:
        return "N"

def ensure_fasta_with_choice(build, signal_callback, ask_callback, stop_event=None):
    """Lädt FASTA bei Bedarf. stop_event prüft auf Abbruch."""
    fasta_path = FASTA_PATHS.get(build)
    if not fasta_path: return None
    fai_path = fasta_path + ".fai"
    
    # Prüfen ob Abbruch
    if stop_event and stop_event.is_set(): return None

    if os.path.exists(fasta_path) and os.path.exists(fai_path):
        return fasta_path
    
    if os.path.exists(fasta_path):
        signal_callback.emit(f"Erstelle Index für {fasta_path}...")
        try:
            build_fasta_index(fasta_path)
            return fasta_path
        except Exception as e:
            signal_callback.emit(f"[ERROR] Indexierung fehlgeschlagen: {e}")
            return None

    # Frage über Signal an GUI
    should_download = ask_callback.emit(
        "Referenz laden?",
        f"Die Referenz-Datei für {build} fehlt (~850 MB).\nHerunterladen? (Empfohlen)"
    )
    
    if not should_download: return None

    url = FASTA_URLS[build]
    gz_path = fasta_path + ".gz"
    try:
        signal_callback.emit(f"Starte Download {build}...")
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(gz_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192*16):
                    if stop_event and stop_event.is_set():
                        signal_callback.emit("Download abgebrochen.")
                        f.close()
                        os.remove(gz_path)
                        return None
                    if chunk:
                        f.write(chunk)
        
        signal_callback.emit("Download fertig. Entpacke Datei...")
        with gzip.open(gz_path, "rb") as f_in, open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)
        
        signal_callback.emit("Erstelle Index...")
        build_fasta_index(fasta_path)
        return fasta_path
    except Exception as e:
        signal_callback.emit(f"Download Fehler: {e}")
        if os.path.exists(gz_path): os.remove(gz_path)
        if os.path.exists(fasta_path): os.remove(fasta_path)
        return None

# -----------------------------
# dbSNP & Parsing
# -----------------------------
def nc_to_chrom(nc):
    ignore = ("XM_", "NM_", "XP_", "NP_", "NG_", "NW_", "NR_", "XR_")
    if nc.startswith(ignore): return None
    if nc.startswith("NC_012920"): return "MT"
    m = re.match(r"NC_(\d{6})\.\d+", nc)
    if m:
        num = int(m.group(1))
        if 1 <= num <= 22: return str(num)
        if num == 23: return "X"
        if num == 24: return "Y"
    return None

def parse_refsnp_payload(rsid, data):
    results = {}
    try: pwa = data["primary_snapshot_data"]["placements_with_allele"]
    except: return results
    for p in pwa:
        alleles = p.get("alleles") or []
        for allele_entry in alleles:
            spdi = allele_entry.get("allele", {}).get("spdi", {})
            seq_id, pos, deleted = spdi.get("seq_id"), spdi.get("position"), spdi.get("deleted_sequence", "").upper()
            if not seq_id or pos is None: continue
            chrom = nc_to_chrom(seq_id)
            if not chrom: continue
            if not deleted: deleted = "-"
            for t in p.get("placement_annot", {}).get("seq_id_traits_by_assembly") or []:
                asm = t.get("assembly_name")
                if asm and (asm.startswith("GRCh37") or asm.startswith("GRCh38")):
                    base = "GRCh37" if asm.startswith("GRCh37") else "GRCh38"
                    results[base] = {"chrom": chrom, "pos": int(pos)+1, "ref": deleted}
    return results

def fetch_single_refsnp(rsid, session):
    rid = rsid.lower().replace("rs", "")
    url = DBSNP_URL + rid
    attempt = 0
    while attempt <= RETRY_MAX:
        try:
            r = session.get(url, timeout=TIMEOUT_SEC)
            if r.status_code == 200:
                try: return rsid, r.json()
                except: return rsid, None
            if r.status_code in (429, 502, 503, 504):
                time.sleep(backoff_time(attempt))
                attempt += 1
                continue
            return rsid, None
        except:
            time.sleep(backoff_time(attempt))
            attempt += 1
    return rsid, None

def adaptive_parallel_fetch(rsids, cache, signal_callback, stop_event):
    to_fetch = [r for r in rsids if r not in cache]
    if not to_fetch: return

    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(pool_connections=MAX_WORKERS, pool_maxsize=MAX_WORKERS)
    session.mount("https://", adapter)
    
    workers = min(MAX_WORKERS, max(MIN_WORKERS, (cpu_count() or 4) * 2))
    processed = 0
    total = len(to_fetch)
    save_cnt = 0

    while to_fetch:
        if stop_event.is_set():
            signal_callback.emit("Abbruch durch Benutzer während Download.")
            break

        chunk_size = max(1, min(len(to_fetch), workers * 6))
        chunk = [to_fetch.pop() for _ in range(chunk_size)]
        t0 = time.time()
        
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(fetch_single_refsnp, r, session): r for r in chunk}
            for fut in as_completed(futures):
                if stop_event.is_set(): break # Schneller Exit
                
                rsid, data = fut.result()
                if data:
                    res = parse_refsnp_payload(rsid, data)
                    for bld, rec in res.items():
                        cache_upsert(cache, rsid, bld, rec["chrom"], rec["pos"], rec["ref"])
                processed += 1
                save_cnt += 1

        if save_cnt >= 500:
            save_cache(cache)
            save_cnt = 0
        
        if not stop_event.is_set():
            dt = time.time() - t0
            speed = len(chunk)/max(0.01, dt)
            signal_callback.emit(f"Download: {processed}/{total} (Speed: {speed:.1f}/s)")
    
    save_cache(cache)
    session.close()

def parse_genotype_file(file_path):
    """Parse DTC genotype files in 23andMe-style TSV format.

    Supports any provider using tab-separated: rsid, chromosome, position, genotype
    Also handles comma-separated (CSV) files from MyHeritage, FTDNA, tellmeGen.
    """
    variants = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith('"') or not line.strip(): continue
            # Try tab first, fall back to comma
            parts = line.strip().split("\t")
            if len(parts) < 4:
                parts = line.strip().split(",")
            if len(parts) < 4: continue
            # Strip quotes from CSV fields
            rsid, chrom, pos, gt = (p.strip().strip('"') for p in parts[:4])
            # Skip header rows
            if rsid.lower() in ("rsid", "snp", "marker", "name"): continue
            chrom = chrom.replace("chr", "").upper()
            if chrom == "M": chrom = "MT"
            try: variants.append((rsid, chrom, int(pos), gt.strip()))
            except: continue
    return variants

def detect_build_robust(variants, cache, signal_callback, stop_event):
    candidates = [v for v in variants if is_rs_id(v[0])][:200]
    adaptive_parallel_fetch([c[0] for c in candidates], cache, signal_callback, stop_event)
    
    if stop_event.is_set(): return None

    m37, m38, tol = 0, 0, 5
    for rsid, _, pos, _ in candidates:
        entry = cache.get(rsid, {}).get("assemblies", {})
        if "GRCh37" in entry and abs(entry["GRCh37"]["pos"] - pos) <= tol: m37 += 1
        if "GRCh38" in entry and abs(entry["GRCh38"]["pos"] - pos) <= tol: m38 += 1
    
    return "GRCh38" if m38 > m37 else "GRCh37"

def create_vcf(variants, build, out_path, cache, fasta_path=None, sex="unknown", signal_callback=None, stop_event=None, progress_signal=None):
    fai_index = load_fai_index(fasta_path + ".fai") if fasta_path else {}
    
    def get_ref(chrom, pos):
        if fasta_path and fai_index:
            b = fetch_base_from_fasta(fasta_path, fai_index, chrom, pos)
            if b and b != "N": return b
        entry = cache.get(rsid, {}).get("assemblies", {}).get(build)
        return entry.get("ref", "N") if entry else "N"

    written = 0
    total_variants = len(variants)
    
    with open(out_path, "w", encoding="utf-8") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##reference={build}\n")
        vcf.write(f"##source=Genotype_to_VCF_Pro_v{APP_VERSION}\n")
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write(f'##INFO=<ID=I_ID,Number=1,Type=String,Description="Original internal ID">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for idx, (rsid, chrom, pos, genotype) in enumerate(variants):
            if idx % 1000 == 0:
                if stop_event and stop_event.is_set():
                    signal_callback.emit("Schreiben abgebrochen.")
                    return 0
                if progress_signal:
                    progress_signal.emit(int((idx / total_variants) * 100))

            genotype = genotype.replace("_", "-")
            if genotype in ("--", "-", "00"): continue 

            ref_base = get_ref(chrom, pos)
            if not ref_base or ref_base in (".", "N", "-"): continue
            ref_base = ref_base.upper()

            ploid = ploidy_for_site(chrom, pos, build, sex)
            if ploid == 0: continue

            alleles = list(genotype)
            FILTER, INFO = "PASS", "."

            if rsid.startswith("i"):
                mapped = lookup_rsid_from_cache(chrom, pos, build, cache)
                if mapped: INFO, rsid = f"I_ID={rsid}", mapped
                else: rsid = "."

            is_simple_snp = all(a in "ACGT" for a in alleles) and ref_base in "ACGT"
            
            if is_simple_snp:
                sample_alts = sorted(list(set([a for a in alleles if a != ref_base])))
                if not sample_alts:
                    ALT, GT = ".", ("0/0" if ploid == 2 else "0")
                else:
                    ALT = ",".join(sample_alts)
                    GT = gt_from_bases_snp(alleles if ploid == 2 else alleles[:1], ref_base, sample_alts, ploid)
                
                if GT:
                    vcf.write(f"{chrom}\t{pos}\t{rsid}\t{ref_base}\t{ALT}\t.\t{FILTER}\t{INFO}\tGT\t{GT}\n")
                    written += 1

    return written

# -----------------------------
# PyQt6 Worker Thread
# -----------------------------
class ConversionWorker(QThread):
    log_signal = pyqtSignal(str)
    progress_signal = pyqtSignal(int)
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(str)
    ask_dialog_signal = pyqtSignal(str, str) # Title, Message -> Returns bool via wait condition?
    
    # Lösung für Dialog: Signal senden, Variable setzen
    dialog_result = None
    dialog_event = threading.Event()

    def __init__(self, file_path, sex, build, cache):
        super().__init__()
        self.file_path = file_path
        self.sex = sex
        self.build_input = build
        self.cache = cache
        self.is_interrupted = threading.Event()

    def stop(self):
        self.is_interrupted.set()

    def ask_wrapper(self, title, msg):
        self.dialog_result = None
        self.dialog_event.clear()
        self.ask_dialog_signal.emit(title, msg)
        self.dialog_event.wait() # Wait for GUI to answer
        return self.dialog_result

    def run(self):
        try:
            self.log_signal.emit("Lade Varianten...")
            variants = parse_genotype_file(self.file_path)
            
            if not variants:
                self.error_signal.emit("Keine gültigen Varianten in Datei gefunden.")
                return

            if self.is_interrupted.is_set(): return

            # Sex Detection
            if self.sex == "Auto":
                self.log_signal.emit("Ermittle Geschlecht (Y-Check)...")
                self.sex = detect_sex_from_variants(variants)
                self.log_signal.emit(f"Geschlecht erkannt: {self.sex}")

            # Build Detection
            if self.build_input == "Auto":
                self.log_signal.emit("Ermittle Build (dbSNP Check)...")
                detected_build = detect_build_robust(variants, self.cache, self.log_signal, self.is_interrupted)
                if not detected_build: 
                    if self.is_interrupted.is_set(): return
                    self.error_signal.emit("Build konnte nicht erkannt werden.")
                    return
                self.build_input = detected_build

            self.log_signal.emit(f"Build: {self.build_input}")
            if self.build_input not in FASTA_PATHS:
                self.error_signal.emit("Ungültiger Build.")
                return

            # FASTA
            fasta_path = ensure_fasta_with_choice(
                self.build_input, self.log_signal, self.ask_wrapper, self.is_interrupted
            )
            
            if self.is_interrupted.is_set(): return
            
            if fasta_path: 
                self.log_signal.emit("Lokale FASTA wird genutzt (Schnellmodus).")
                self.log_signal.emit("Überspringe dbSNP-Download für SNPs.")
            else: 
                self.log_signal.emit("Keine FASTA. Fallback Modus (langsamer).")
                # Cache Fill nur wenn keine FASTA da
                missing = [v[0] for v in variants if is_rs_id(v[0]) and v[0] not in self.cache]
                if missing:
                    self.log_signal.emit(f"Lade {len(missing)} Metadaten nach...")
                    adaptive_parallel_fetch(missing, self.cache, self.log_signal, self.is_interrupted)

            if self.is_interrupted.is_set(): return

            # Convert
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            out_name = f"{os.path.splitext(self.file_path)[0]}_{self.build_input}_{ts}.vcf"
            
            self.log_signal.emit(f"Schreibe VCF: {os.path.basename(out_name)}...")
            count = create_vcf(
                variants, self.build_input, out_name, self.cache, 
                fasta_path, self.sex, self.log_signal, self.is_interrupted, self.progress_signal
            )
            
            if not self.is_interrupted.is_set():
                self.progress_signal.emit(100)
                self.finished_signal.emit(f"Erfolg! {count} Varianten geschrieben.\nDatei: {out_name}")
            else:
                self.log_signal.emit("Vorgang abgebrochen.")

        except Exception as e:
            self.error_signal.emit(str(e))

# -----------------------------
# Modern GUI Class
# -----------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_TITLE)
        self.resize(900, 650)
        self.cache = load_cache()
        self.file_path = None
        self.worker = None
        
        self.init_ui()
        self.apply_styles()

    def init_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)

        # --- Header Section ---
        header_layout = QHBoxLayout()
        
        self.btn_file = QPushButton("📁 Genotype-Datei öffnen")
        self.btn_file.setMinimumHeight(40)
        self.btn_file.clicked.connect(self.select_file)
        
        self.lbl_file = QLabel("Keine Datei ausgewählt")
        self.lbl_file.setStyleSheet("color: #888; font-style: italic;")
        
        header_layout.addWidget(self.btn_file)
        header_layout.addWidget(self.lbl_file, 1)
        layout.addLayout(header_layout)

        # --- Options Section ---
        options_frame = QFrame()
        options_frame.setObjectName("OptionsFrame")
        options_layout = QHBoxLayout(options_frame)
        
        # Sex
        options_layout.addWidget(QLabel("Geschlecht:"))
        self.combo_sex = QComboBox()
        self.combo_sex.addItems(["Auto", "female", "male"])
        options_layout.addWidget(self.combo_sex)
        
        # Build
        options_layout.addWidget(QLabel("Build:"))
        self.combo_build = QComboBox()
        self.combo_build.addItems(["Auto", "GRCh37", "GRCh38"])
        options_layout.addWidget(self.combo_build)
        
        options_layout.addStretch()
        layout.addWidget(options_frame)

        # --- Action Buttons ---
        btn_layout = QHBoxLayout()
        
        self.btn_start = QPushButton("🚀 Konvertierung Starten")
        self.btn_start.setMinimumHeight(50)
        self.btn_start.setObjectName("StartButton")
        self.btn_start.clicked.connect(self.start_process)
        
        self.btn_cancel = QPushButton("🛑 Abbrechen")
        self.btn_cancel.setMinimumHeight(50)
        self.btn_cancel.setObjectName("CancelButton")
        self.btn_cancel.setEnabled(False)
        self.btn_cancel.clicked.connect(self.cancel_process)

        btn_layout.addWidget(self.btn_start, 2)
        btn_layout.addWidget(self.btn_cancel, 1)
        layout.addLayout(btn_layout)

        # --- Progress ---
        self.progress_bar = QProgressBar()
        self.progress_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.progress_bar)

        # --- Log ---
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        self.txt_log.setFont(QFont("Consolas", 9))
        layout.addWidget(self.txt_log)

    def apply_styles(self):
        # Fusion Style Base
        QApplication.setStyle(QStyleFactory.create("Fusion"))
        
        # Dark Palette
        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
        palette.setColor(QPalette.ColorRole.WindowText, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.Base, QColor(25, 25, 25))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
        palette.setColor(QPalette.ColorRole.ToolTipBase, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.ToolTipText, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
        palette.setColor(QPalette.ColorRole.ButtonText, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.BrightText, Qt.GlobalColor.red)
        palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
        palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
        palette.setColor(QPalette.ColorRole.HighlightedText, Qt.GlobalColor.black)
        QApplication.setPalette(palette)

        # Stylesheet for specific widgets
        self.setStyleSheet("""
            QPushButton {
                background-color: #444;
                border: 1px solid #555;
                border-radius: 4px;
                padding: 5px;
                color: white;
            }
            QPushButton:hover {
                background-color: #555;
            }
            QPushButton:pressed {
                background-color: #333;
            }
            QPushButton#StartButton {
                background-color: #2e7d32; /* Green */
                font-weight: bold;
                font-size: 14px;
            }
            QPushButton#StartButton:hover { background-color: #388e3c; }
            
            QPushButton#CancelButton {
                background-color: #c62828; /* Red */
                font-weight: bold;
            }
            QPushButton#CancelButton:hover { background-color: #d32f2f; }
            QPushButton:disabled {
                background-color: #333;
                color: #777;
            }

            QFrame#OptionsFrame {
                background-color: #444;
                border-radius: 5px;
            }
            QProgressBar {
                border: 2px solid grey;
                border-radius: 5px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #2e7d32;
                width: 20px;
            }
        """)

    def select_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Genotype-Datei wählen", "", "Genotype Files (*.txt *.csv *.tsv);;All Files (*)")
        if fname:
            self.file_path = fname
            self.lbl_file.setText(os.path.basename(fname))
            self.lbl_file.setStyleSheet("color: white; font-weight: bold;")

    def log(self, msg):
        self.txt_log.append(msg)

    def handle_ask_dialog(self, title, msg):
        reply = QMessageBox.question(self, title, msg, 
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, 
                                     QMessageBox.StandardButton.Yes)
        # Ergebnis an Worker zurückgeben
        if self.worker:
            self.worker.dialog_result = (reply == QMessageBox.StandardButton.Yes)
            self.worker.dialog_event.set()

    def start_process(self):
        if not self.file_path:
            QMessageBox.warning(self, "Fehler", "Bitte zuerst eine Datei auswählen.")
            return

        self.txt_log.clear()
        self.progress_bar.setValue(0)
        self.btn_start.setEnabled(False)
        self.btn_file.setEnabled(False)
        self.btn_cancel.setEnabled(True)
        
        # Worker Setup
        sex = self.combo_sex.currentText()
        build = self.combo_build.currentText()
        
        self.worker = ConversionWorker(self.file_path, sex, build, self.cache)
        self.worker.log_signal.connect(self.log)
        self.worker.progress_signal.connect(self.progress_bar.setValue)
        self.worker.error_signal.connect(lambda e: QMessageBox.critical(self, "Fehler", e))
        self.worker.finished_signal.connect(lambda m: QMessageBox.information(self, "Fertig", m))
        
        # Cleanup nach Ende (egal ob Erfolg oder Fehler)
        self.worker.finished.connect(self.process_finished)
        self.worker.ask_dialog_signal.connect(self.handle_ask_dialog)
        
        self.worker.start()

    def cancel_process(self):
        if self.worker and self.worker.isRunning():
            self.log("🛑 Abbruch angefordert...")
            self.worker.stop()
            # Wir lassen den Thread 'gracefully' auslaufen, daher kein terminate()
            self.btn_cancel.setEnabled(False)

    def process_finished(self):
        self.btn_start.setEnabled(True)
        self.btn_file.setEnabled(True)
        self.btn_cancel.setEnabled(False)
        self.worker = None

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())