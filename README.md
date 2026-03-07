# Genotype-to-VCF Pro Converter

Eine Desktop-Anwendung zur Konvertierung von DTC (Direct-to-Consumer) DNA-Rohdaten in das standardisierte **VCF 4.2** Format. Mit moderner PyQt6-Oberflaeche, Unterstuetzung fuer GRCh37 und GRCh38 Referenzgenome und automatischer Build-Erkennung.

Urspruenglich fuer 23andMe-Exporte entwickelt, funktioniert es mit **jedem Anbieter**, der das gleiche Tab-separierte Format verwendet (`rsid  chromosome  position  genotype`).

## Funktionen

| Funktion | Beschreibung |
|---|---|
| **Dual-Referenzgenom** | GRCh37 (hg19) und GRCh38 (hg38) |
| **Auto Build-Erkennung** | Erkennt Genom-Build via dbSNP-Positionsvalidierung |
| **Auto Geschlechtserkennung** | Bestimmt biologisches Geschlecht anhand Y-Chromosom-Varianten |
| **PAR-Region-Behandlung** | Korrekte Ploidie fuer pseudoautosomale Regionen auf X/Y |
| **dbSNP-Integration** | NCBI REST API fuer rsID-Abfragen und REF-Basen |
| **Persistenter Cache** | Lokaler Cache fuer schnelle wiederholte Konvertierungen |
| **Adaptives Threading** | 4-200 Worker-Threads, Ziel 70% CPU-Auslastung |
| **FASTA-Referenz** | Optionale lokale FASTA-Datei fuer Offline-REF-Abfrage |
| **Moderne GUI** | PyQt6 Dark Theme mit Fortschrittsanzeige und Abbruch-Option |

## Unterstuetzte Eingabeformate

Das Tool liest Tab-separierte Dateien (TSV) mit vier Spalten:

```
# rsid  chromosome  position  genotype
rs12564807	1	734462	AA
rs3131972	1	752721	AG
```

Zeilen, die mit `#` beginnen, werden als Kommentare uebersprungen.

### Getestete Anbieter

| Anbieter | Kompatibel | Hinweise |
|---|---|---|
| **23andMe** (v3/v4/v5) | Ja | Natives Format |
| **Genes for Good** | Ja | Exportiert im 23andMe-Format |
| **Mapmygenome** | Ja | 23andMe-kompatibles Format |
| **MyHeritage** | Ja | CSV mit 4 Spalten (automatisch erkannt) |
| **Family Tree DNA** | Ja | CSV mit 4 Spalten (automatisch erkannt) |
| **tellmeGen** | Ja | CSV mit 4 Spalten (automatisch erkannt) |
| **AncestryDNA** | Nein | Verwendet 5 Spalten (allele1, allele2 getrennt) |
| **LivingDNA** | Nein | Andere Spaltenreihenfolge |

> **Tipp:** TSV (Tab-separiert) und CSV (Komma-separiert) werden automatisch erkannt. Jede Datei mit vier Spalten (`rsid, chrom, pos, genotype`) funktioniert, unabhaengig vom Anbieter.

## Installation

### Option 1: Windows Executable (kein Python noetig)

Die aktuelle `23toVCF_Pro.exe` von der [Releases](../../releases)-Seite herunterladen und direkt starten.

### Option 2: Aus dem Quellcode

**Voraussetzungen:** Python 3.8+

```bash
git clone https://github.com/lukisch/genotype-to-vcf.git
cd genotype-to-vcf
pip install -r requirements.txt
python Make23toVCF3.py
```

### Option 3: Eigene EXE erstellen

```bash
pip install pyinstaller
pyinstaller --onefile --windowed --name 23toVCF_Pro Make23toVCF3.py
```

## Verwendung

1. Anwendung starten
2. **"Open File"** klicken und Rohdatendatei (`.txt`) auswaehlen
3. **Geschlecht** (`Auto` / `female` / `male`) und **Build** (`Auto` / `GRCh37` / `GRCh38`) waehlen
4. **"Start Conversion"** klicken
5. Die VCF-Datei wird neben der Eingabedatei gespeichert

### Erster Start

Beim ersten Start ohne lokale FASTA-Referenz:
- Verwendung der **NCBI dbSNP API** fuer Referenzbasen-Abfrage (langsamer, Internet erforderlich)
- Angebot zum **FASTA-Download** (~850 MB pro Build) fuer schnellere Offline-Konvertierungen
- Aufbau eines **lokalen Caches** (`cache.json`) fuer alle folgenden Konvertierungen

### Konvertierungs-Pipeline

```
Eingabedatei (.txt)
    |
    v
TSV parsen (rsid, chrom, pos, genotype)
    |
    v
Build erkennen (GRCh37 vs GRCh38) via dbSNP-Validierung
    |
    v
Geschlecht erkennen (Y-Chromosom-Varianten-Anzahl)
    |
    v
REF-Basen aufloesen (FASTA > Cache > dbSNP API > ueberspringen)
    |
    v
VCF 4.2 schreiben mit korrekter Ploidie und Genotyp-Calls
    |
    v
Ausgabe: sample_GRCh37_20260213_143000.vcf
```

## VCF-Ausgabeformat

```vcf
##fileformat=VCFv4.2
##reference=GRCh37
##source=23andMe_Pro_Converter
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=I_ID,Number=1,Type=String,Description="Original internal ID">
#CHROM  POS     ID          REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
1       734462  rs12564807  A    .    .     PASS    .     GT      0/0
1       752721  rs3131972   A    G    .     PASS    .     GT      0/1
```

### Genotyp-Kodierung

| Ploidie | Kontext | Beispiel |
|---|---|---|
| Diploid (0/0, 0/1, 1/1) | Autosomen, weibliches X, PAR-Regionen | `0/1` |
| Haploid (0, 1) | Maennliches X (non-PAR), maennliches Y (non-PAR), MT | `1` |
| Uebersprungen | Weibliche Y-Varianten | - |

## Funktionsweise

### Build-Erkennung

Das Tool nimmt bis zu 200 Varianten mit rsIDs und fragt die NCBI dbSNP API nach deren genomischen Positionen auf GRCh37 und GRCh38. Der Build mit den meisten Positionsuebereinstimmungen (Toleranz 5 bp) wird gewaehlt.

### REF-Basen-Aufloesung

Referenzbasen werden in dieser Prioritaet aufgeloest:

1. **Lokale FASTA** - Byte-exakter Lookup via `.fai`-Index (schnellste)
2. **Lokaler Cache** - Zuvor abgerufene dbSNP-Daten
3. **dbSNP API** - Live NCBI REST API Abfrage
4. **Ueberspringen** - Varianten ohne aufgeloeste REF-Base werden ausgeschlossen

### Caching

Der persistente `cache.json` speichert dbSNP API-Antworten mit Zeitstempeln. Folgekonvertierungen von Dateien mit ueberlappenden SNPs sind deutlich schneller. Der Cache nutzt atomare Schreibvorgaenge mit File-Locking fuer Thread-Sicherheit.

## Technische Details

- **Sprache:** Python 3.8+
- **GUI:** PyQt6 mit Fusion Dark Theme
- **Bioinformatik:** pyfaidx fuer FASTA-Indexierung
- **API:** NCBI dbSNP REST API (`https://api.ncbi.nlm.nih.gov/variation/v0/`)
- **Threading:** `ThreadPoolExecutor` mit CPU-adaptiver Worker-Anzahl
- **VCF-Standard:** v4.2 ([Spezifikation](https://samtools.github.io/hts-specs/VCFv4.2.pdf))

## Datenschutz

Dieses Tool verarbeitet genetische Daten lokal auf Ihrem Rechner. Keine Daten werden an externe Server gesendet, ausser:
- **NCBI dbSNP API** Abfragen, die ausschliesslich rsIDs enthalten (z.B. `rs12345`) zur Positionsaufloesung
- **Ensembl FTP** fuer optionale FASTA-Referenzgenom-Downloads

Genotyp-Daten, persoenliche Identifikatoren oder Rohdateien werden niemals uebertragen.

---

## English

A desktop application for converting DTC (Direct-to-Consumer) DNA raw data files into the standardized **VCF 4.2** format. Built with a modern PyQt6 GUI, it supports both GRCh37 and GRCh38 reference genomes with automatic build detection.

Originally designed for 23andMe exports, it works with **any provider** that uses the same tab-separated format (`rsid  chromosome  position  genotype`).

### Features

| Feature | Description |
|---|---|
| **Dual Reference Genome** | GRCh37 (hg19) and GRCh38 (hg38) |
| **Auto Build Detection** | Detects genome build via dbSNP position validation |
| **Auto Sex Detection** | Infers biological sex from Y chromosome variants |
| **PAR Region Handling** | Correct ploidy for pseudo-autosomal regions on X/Y |
| **dbSNP Integration** | NCBI REST API for rsID lookup and REF base retrieval |
| **Persistent Cache** | Local cache for fast repeated conversions |
| **Adaptive Threading** | 4-200 worker threads, targeting 70% CPU usage |
| **FASTA Reference** | Optional local FASTA for offline REF base lookup |
| **Modern GUI** | PyQt6 dark theme with progress tracking and cancel support |

### Supported Input Formats

This tool reads tab-separated files (TSV) with four columns:

```
# rsid  chromosome  position  genotype
rs12564807	1	734462	AA
rs3131972	1	752721	AG
```

Lines starting with `#` are treated as comments and skipped.

### Tested Providers

| Provider | Compatible | Notes |
|---|---|---|
| **23andMe** (v3/v4/v5) | Yes | Native format |
| **Genes for Good** | Yes | Exports in 23andMe format |
| **Mapmygenome** | Yes | Uses 23andMe-compatible format |
| **MyHeritage** | Yes | CSV with 4 columns (auto-detected) |
| **Family Tree DNA** | Yes | CSV with 4 columns (auto-detected) |
| **tellmeGen** | Yes | CSV with 4 columns (auto-detected) |
| **AncestryDNA** | No | Uses 5 columns (allele1, allele2 separate) |
| **LivingDNA** | No | Different column order |

> **Tip:** Both TSV (tab-separated) and CSV (comma-separated) files are auto-detected. Any file with four columns (`rsid, chrom, pos, genotype`) will work, regardless of the provider.

### Installation

#### Option 1: Windows Executable (No Python Required)

Download the latest `23toVCF_Pro.exe` from the [Releases](../../releases) page and run it directly.

#### Option 2: From Source

**Requirements:** Python 3.8+

```bash
git clone https://github.com/lukisch/genotype-to-vcf.git
cd genotype-to-vcf
pip install -r requirements.txt
python Make23toVCF3.py
```

#### Option 3: Build Your Own Executable

```bash
pip install pyinstaller
pyinstaller --onefile --windowed --name 23toVCF_Pro Make23toVCF3.py
```

### Usage

1. Launch the application
2. Click **"Open File"** and select your raw data file (`.txt`)
3. Choose **Sex** (`Auto` / `female` / `male`) and **Build** (`Auto` / `GRCh37` / `GRCh38`)
4. Click **"Start Conversion"**
5. The VCF file is saved alongside the input file

#### First Run

On first run without a local FASTA reference, the tool will:
- Use the **NCBI dbSNP API** to look up reference bases (slower, requires internet)
- Offer to **download the FASTA reference** (~850 MB per build) for faster offline conversions
- Build a **local cache** (`cache.json`) that speeds up all subsequent conversions

#### Conversion Pipeline

```
Input File (.txt)
    |
    v
Parse TSV (rsid, chrom, pos, genotype)
    |
    v
Detect Build (GRCh37 vs GRCh38) via dbSNP validation
    |
    v
Detect Sex (Y chromosome variant count)
    |
    v
Resolve REF bases (FASTA > Cache > dbSNP API > skip)
    |
    v
Write VCF 4.2 with correct ploidy and genotype calls
    |
    v
Output: sample_GRCh37_20260213_143000.vcf
```

### VCF Output Format

```vcf
##fileformat=VCFv4.2
##reference=GRCh37
##source=23andMe_Pro_Converter
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=I_ID,Number=1,Type=String,Description="Original internal ID">
#CHROM  POS     ID          REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
1       734462  rs12564807  A    .    .     PASS    .     GT      0/0
1       752721  rs3131972   A    G    .     PASS    .     GT      0/1
```

#### Genotype Encoding

| Ploidy | Context | Example |
|---|---|---|
| Diploid (0/0, 0/1, 1/1) | Autosomes, female X, PAR regions | `0/1` |
| Haploid (0, 1) | Male X (non-PAR), male Y (non-PAR), MT | `1` |
| Skipped | Female Y variants | - |

### How It Works

#### Build Detection

The tool samples up to 200 variants with rsIDs and queries the NCBI dbSNP API for their genomic positions on both GRCh37 and GRCh38. The build with the most position matches (within 5 bp tolerance) is selected.

#### REF Base Resolution

Reference bases are resolved in this priority order:

1. **Local FASTA** - Byte-exact lookup via `.fai` index (fastest)
2. **Local Cache** - Previously fetched dbSNP data
3. **dbSNP API** - Live NCBI REST API query
4. **Skip** - Variants without a resolved REF base are excluded

#### Caching

The persistent `cache.json` stores dbSNP API responses with timestamps. Subsequent conversions of files with overlapping SNPs are significantly faster. The cache uses atomic writes with file locking for thread safety.

### Technical Details

- **Language:** Python 3.8+
- **GUI:** PyQt6 with Fusion dark theme
- **Bioinformatics:** pyfaidx for FASTA indexing
- **API:** NCBI dbSNP REST API (`https://api.ncbi.nlm.nih.gov/variation/v0/`)
- **Threading:** `ThreadPoolExecutor` with CPU-adaptive worker count
- **VCF Standard:** v4.2 ([specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf))

### Privacy

This tool processes genetic data locally on your machine. No data is sent to external servers except:
- **NCBI dbSNP API** queries containing only rsIDs (e.g., `rs12345`) to resolve reference positions
- **Ensembl FTP** for optional FASTA reference genome downloads

No genotype data, personal identifiers, or raw files are ever transmitted.

### Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Acknowledgments

- [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/) for the variant database API
- [Ensembl](https://www.ensembl.org/) for reference genome sequences
- [pyfaidx](https://github.com/mdshw5/pyfaidx) for FASTA indexing

## License

[MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.
