# Minimal AlphaSimR Functions and Classes Required for Each Tutorial

This document lists the minimal set of AlphaSimR functions and classes needed to run each tutorial in this repository.

## 01_LineBreeding

### 01_PhenotypicSelection/01_MassSelection
**Functions:**
- `runMacs()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `setPheno()` - Set phenotypes
- `selectInd()` - Select individuals
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation (base R, used with AlphaSimR objects)

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 01_PhenotypicSelection/02_SingleSeedDescent
**Functions:**
- `runMacs()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `self()` - Self-pollination
- `setPheno()` - Set phenotypes
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 01_PhenotypicSelection/03_PedigreeSelection
**Functions:**
- `runMacs()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `self()` - Self-pollination
- `setPheno()` - Set phenotypes
- `selectInd()` - Select individuals
- `meanP()` - Calculate mean phenotype
- `mergePops()` - Merge populations
- `nInd()` - Get number of individuals
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `accuracy_family()` - Custom function (uses `meanP()`, `meanG()`, `cor()`)

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 01_PhenotypicSelection/04_DoubledHaploid
**Functions:**
- `runMacs()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `makeDH()` - Create doubled haploids
- `setPheno()` - Set phenotypes
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 02_GenomicSelection
**Functions:**
- `runMacs()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `makeDH()` - Create doubled haploids
- `setPheno()` - Set phenotypes
- `setEBV()` - Set estimated breeding values
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `RRBLUP()` - Ridge regression BLUP model
- `c()` - Combine populations (base R, used with AlphaSimR objects)
- `@nInd` - Access number of individuals slot

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 03_TwoPartGS
**Functions:**
- `runMacs()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `makeDH()` - Create doubled haploids
- `setPheno()` - Set phenotypes
- `setEBV()` - Set estimated breeding values
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `RRBLUP()` - Ridge regression BLUP model
- `c()` - Combine populations
- `@nInd` - Access number of individuals slot

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

## 02_ClonalBreeding

### 01_PhenotypicSelection
**Functions:**
- `runMacs2()` - Generate initial haplotypes (clonal species)
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `setPheno()` - Set phenotypes (with `p` parameter for year effects)
- `selectInd()` - Select individuals
- `gv()` - Get genetic values
- `pheno()` - Get phenotypes
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 02_PedigreeSelection
**Functions:**
- `runMacs2()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `setPheno()` - Set phenotypes (with `p` parameter)
- `setEBV()` - Set estimated breeding values (from pedigree model)
- `selectInd()` - Select individuals
- `gv()` - Get genetic values
- `pheno()` - Get phenotypes
- `ebv()` - Get estimated breeding values
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `nInd()` - Get number of individuals
- `@fixEff` - Access fixed effects slot
- `@id` - Access individual IDs slot
- `@father` - Access father IDs slot
- `@mother` - Access mother IDs slot
- `@nInd` - Access number of individuals slot

**Note:** This tutorial also uses external packages (asreml) for pedigree BLUP, but the minimal AlphaSimR requirement is `setEBV()` which can work with internal AlphaSimR solver.

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 03_GenomicSelection
**Functions:**
- `runMacs2()` - Generate initial haplotypes
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE
- `SP$setTrackPed()` - Enable pedigree tracking
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `setPheno()` - Set phenotypes (with `p` parameter)
- `setEBV()` - Set estimated breeding values (from genomic model)
- `selectInd()` - Select individuals
- `gv()` - Get genetic values
- `pheno()` - Get phenotypes
- `ebv()` - Get estimated breeding values
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `RRBLUP()` - Ridge regression BLUP model
- `nInd()` - Get number of individuals
- `@fixEff` - Access fixed effects slot
- `@nInd` - Access number of individuals slot

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

## 03_HybridBreeding

### 01_PhenotypicSelection
**Functions:**
- `runMacs()` - Generate initial haplotypes (with `split` parameter for heterotic pools)
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE
- `SP$setTrackPed()` - Enable pedigree tracking
- `SP$setVarE()` - Set permanent error variance
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `makeDH()` - Create doubled haploids
- `setPheno()` - Set phenotypes (with `p` parameter)
- `setPhenoGCA()` - Set phenotypes using general combining ability
- `selectInd()` - Select individuals
- `hybridCross()` - Create hybrid crosses
- `calcGCA()` - Calculate general combining ability
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `@id` - Access individual IDs slot
- `@mother` - Access mother IDs slot
- `@pheno` - Access phenotypes slot
- `@gv` - Access genetic values slot

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 02_GenomicSelection
**Functions:**
- `runMacs()` - Generate initial haplotypes (with `split` parameter)
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE
- `SP$setTrackPed()` - Enable pedigree tracking
- `SP$setVarE()` - Set permanent error variance
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `makeDH()` - Create doubled haploids
- `setPheno()` - Set phenotypes (with `p` parameter)
- `setPhenoGCA()` - Set phenotypes using general combining ability
- `setEBV()` - Set estimated breeding values
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families
- `hybridCross()` - Create hybrid crosses
- `calcGCA()` - Calculate general combining ability
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `RRBLUP()` - Ridge regression BLUP model
- `@id` - Access individual IDs slot
- `@mother` - Access mother IDs slot
- `@pheno` - Access phenotypes slot
- `@gv` - Access genetic values slot
- `@ebv` - Access estimated breeding values slot

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

### 03_TwoPartGS
**Functions:**
- `runMacs()` - Generate initial haplotypes (with `split` parameter)
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip (if nSnp > 0)
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE
- `SP$setTrackPed()` - Enable pedigree tracking
- `SP$setVarE()` - Set permanent error variance
- `newPop()` - Create founder population
- `randCross()` - Make random crosses
- `makeDH()` - Create doubled haploids
- `setPheno()` - Set phenotypes (with `p` parameter)
- `setPhenoGCA()` - Set phenotypes using general combining ability
- `setEBV()` - Set estimated breeding values
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families
- `hybridCross()` - Create hybrid crosses
- `calcGCA()` - Calculate general combining ability
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation
- `RRBLUP()` - Ridge regression BLUP model
- `@id` - Access individual IDs slot
- `@mother` - Access mother IDs slot
- `@pheno` - Access phenotypes slot
- `@gv` - Access genetic values slot
- `@ebv` - Access estimated breeding values slot

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit)

---

## Summary of Core Functions Across All Tutorials

### Population Creation
- `runMacs()` - For inbred species (wheat, maize)
- `runMacs2()` - For outcrossing/clonal species (tea)
- `newPop()` - Create population from founder haplotypes

### Simulation Parameters
- `SimParam$new()` - Initialize simulation parameters
- `SP$restrSegSites()` - Restrict segregating sites
- `SP$addSnpChip()` - Add SNP chip
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$addTraitADG()` - Add additive-dominance genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking
- `SP$setVarE()` - Set error variance (hybrid breeding)

### Mating
- `randCross()` - Random crosses
- `self()` - Self-pollination
- `makeDH()` - Create doubled haploids
- `hybridCross()` - Create hybrid crosses

### Phenotyping and Selection
- `setPheno()` - Set phenotypes
- `setPhenoGCA()` - Set phenotypes using GCA (hybrid breeding)
- `selectInd()` - Select individuals
- `selectWithinFam()` - Select within families

### Genomic Selection
- `setEBV()` - Set estimated breeding values
- `RRBLUP()` - Ridge regression BLUP model

### Analysis
- `meanG()` - Mean genetic value
- `varG()` - Genetic variance
- `meanP()` - Mean phenotype
- `gv()` - Get genetic values
- `pheno()` - Get phenotypes
- `ebv()` - Get estimated breeding values
- `calcGCA()` - Calculate general combining ability
- `nInd()` - Number of individuals
- `mergePops()` - Merge populations
- `cor()` - Correlation (base R, used with AlphaSimR objects)

### Population Access (Slots)
- `@id` - Individual IDs
- `@mother` - Mother IDs
- `@father` - Father IDs
- `@gv` - Genetic values
- `@pheno` - Phenotypes
- `@ebv` - Estimated breeding values
- `@nInd` - Number of individuals
- `@fixEff` - Fixed effects

### Classes
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit, returned by most functions)

