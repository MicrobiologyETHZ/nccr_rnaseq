# NCCR RNA-seq Pipeline Development Progress

## June 19, 2025 - Functional Annotation Integration

**Goal**: Add comprehensive functional annotation capabilities (eggNOG, dbCAN, cayman) to the NCCR RNA-seq pipeline

### What We Accomplished:

#### 1. Created New Functional Annotation Module
- **File**: `/nccr_rnaseq/rules/functional_annotation.smk`
- **Tools integrated**:
  - **eggNOG-mapper**: Protein functional annotation using Diamond search against eggNOG database
  - **dbCAN**: Carbohydrate-active enzyme (CAZyme) annotation for metabolic pathway analysis
  - **cayman**: Additional functional annotation using CAYMAN v3 database

#### 2. Pipeline Integration
- Added `include: "rules/functional_annotation.smk"` to main Snakefile
- Created master `functional_annotation` rule that coordinates all three annotation tools
- Integrated file dependency tracking with `.done` files for robust workflow management

#### 3. CLI Enhancement
- **New command**: `nccrRna annotate -c config.yaml`
- Supports all standard pipeline options (--local, --dry, --jobs, --partition)
- Runs all three annotation tools in parallel for efficient processing

#### 4. Database Configuration
- **eggNOG**: `/science/ansintsova/eggnog-data/`
- **dbCAN**: `/science/ansintsova/dbCAN/db`
- **cayman**: `/nfs/cds-1/IMB-databases/CAYMAN/v3`
- Added database paths to `configs/rnaseq_config.yaml`

#### 5. Technical Implementation Details
**Command structures implemented**:
```bash
# eggNOG
emapper.py -o {sample} -m diamond --data_dir {db} -i {input.faa} --cpu 16

# dbCAN  
run_dbcan {input.faa} protein -c {input.gff} --dbcan_thread 16 --cgc_substrate

# cayman
cayman annotate_proteome --cutoffs v3/cutoffs.csv v3 {input.faa} --threads 16
```

**Environment**: All tools configured to use `funanot` conda environment

#### 6. Input Requirements
- **Protein sequences**: `{sample}.faa` format
- **Genome annotations**: `{sample}.gff3` format
- **Sample sheet**: CSV file mapping sample names to file paths

### Current Status: ✅ Ready for Testing

### Files Modified:
- `nccr_rnaseq/rules/functional_annotation.smk` (new)
- `nccr_rnaseq/Snakefile` (integration)
- `nccr_rnaseq/main.py` (CLI command)
- `configs/rnaseq_config.yaml` (database paths)

### Usage Example:
```bash
# Generate sample sheet for functional annotation
# Input files: Ass_004_S2.faa, Ass_004_S2.gff3

# Run complete functional annotation
nccrRna annotate -c configs/rnaseq_config.yaml --jobs 5
```

---