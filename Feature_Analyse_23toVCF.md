# Feature-Analyse: 23andMe to VCF Pro Converter

## Kurzbeschreibung
Ein Bioinformatik-Tool zur Konvertierung von 23andMe Genotypisierungsdaten in das standardisierte VCF-Format (Variant Call Format). Unterstützt beide Referenzgenome (GRCh37/GRCh38) und verwendet dbSNP für rsID-Mapping.

---

## ✨ Highlights

| Feature | Beschreibung |
|---------|-------------|
| **Dual-Genome** | GRCh37 (hg19) und GRCh38 (hg38) |
| **Auto-Download** | Referenz-FASTA automatisch laden |
| **dbSNP-Lookup** | NCBI API für rsID-Positionen |
| **Cache-System** | Lokaler Cache für schnelle Lookups |
| **Multi-Threading** | CPU-adaptive Parallelisierung |
| **Geschlechts-Erkennung** | Aus Y-Chromosom-Varianten |
| **Build-Erkennung** | Automatisch aus Datei-Header |
| **PyQt6 GUI** | Moderne Benutzeroberfläche |

---

## 📊 Feature-Vergleich

| Feature | 23toVCF | plink | BCFtools | online Tools |
|---------|:-------:|:-----:|:--------:|:------------:|
| 23andMe-spezifisch | ✅ | ⚠️ | ❌ | ✅ |
| Dual-Genome | ✅ | ✅ | ✅ | ⚠️ |
| dbSNP-Lookup | ✅ | ❌ | ⚠️ | ⚠️ |
| GUI | ✅ | ❌ | ❌ | ✅ |
| Offline-fähig | ✅ | ✅ | ✅ | ❌ |
| Cache | ✅ | N/A | N/A | N/A |
| Auto-Threading | ✅ | ⚠️ | ⚠️ | N/A |

---

## 🎯 Bewertung

### Aktueller Stand: **Production Ready (80%)**

| Kategorie | Bewertung |
|-----------|:---------:|
| Funktionsumfang | ⭐⭐⭐⭐ |
| Wissenschaftlichkeit | ⭐⭐⭐⭐⭐ |
| Performance | ⭐⭐⭐⭐ |

**Gesamtbewertung: 8/10** - Spezialisiertes Bioinformatik-Tool

---

## 🚀 Empfohlene Erweiterungen

1. **Ancestry-Support** - AncestryDNA Format
2. **VCF-Merge** - Mehrere Dateien kombinieren
3. **Annotation** - ClinVar/dbSNP Annotation
4. **Quality-Scores** - Genotyp-Qualität berechnen

---

## 💻 Technische Details

```
Framework:      PyQt6
Bioinformatik:  pyfaidx (FASTA-Index)
API:            NCBI dbSNP REST API
Threading:      ThreadPoolExecutor
CPU-Nutzung:    70% Target, 4-200 Worker
Dateigröße:     775 Zeilen Python
```

### Referenz-Genome:
- GRCh37: ~3 GB FASTA
- GRCh38: ~3 GB FASTA

---

## 🔬 Wissenschaftlicher Kontext

Das VCF-Format ist der Standard für genetische Varianten und wird benötigt für:
- Imputation (TopMed, Michigan Server)
- Ancestry-Analyse (Admixture)
- Krankheitsrisiko-Berechnung
- Forschungsdatenbanken

---
*Analyse erstellt: 02.01.2026*
