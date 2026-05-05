# AlphaSimR Tutorial Roadmap: Progressive Function Introduction

This roadmap shows the progressive introduction of AlphaSimR functions across tutorials, starting with the simplest tutorial. Each tutorial section lists only the **new functions** that haven't been used in previous tutorials.

---

## Tutorial 1: 01_LineBreeding/01_PhenotypicSelection/01_MassSelection

**Foundation Functions (Starting Point):**

**Population Creation:**
- `runMacs()` - Generate initial haplotypes for inbred species
- `newPop()` - Create founder population from haplotypes

**Simulation Parameters:**
- `SimParam$new()` - Create simulation parameters object
- `SP$restrSegSites()` - Restrict segregating sites (QTL vs SNP)
- `SP$addSnpChip()` - Add SNP chip for genotyping
- `SP$addTraitAG()` - Add additive genetic trait
- `SP$setTrackPed()` - Enable pedigree tracking

**Mating:**
- `randCross()` - Make random crosses between parents

**Phenotyping:**
- `setPheno()` - Set phenotypes with error variance

**Selection:**
- `selectInd()` - Select individuals based on criteria

**Analysis:**
- `meanG()` - Calculate mean genetic value
- `varG()` - Calculate genetic variance
- `cor()` - Calculate correlation (base R, used with AlphaSimR objects)

**Classes:**
- `SimParam` - Simulation parameters class
- `Pop` - Population class (implicit, returned by functions)

---

## Tutorial 2: 01_LineBreeding/01_PhenotypicSelection/02_SingleSeedDescent

**New Functions Introduced:**

**Mating:**
- `self()` - Self-pollination to advance generations

**Selection:**
- `selectWithinFam()` - Select individuals within families

---

## Tutorial 3: 01_LineBreeding/01_PhenotypicSelection/03_PedigreeSelection

**New Functions Introduced:**

**Analysis:**
- `meanP()` - Calculate mean phenotype
- `mergePops()` - Merge multiple populations into one
- `nInd()` - Get number of individuals in a population

**Custom Functions:**
- `accuracy_family()` - Custom function to calculate between-family selection accuracy (uses `meanP()`, `meanG()`, `cor()`)

---

## Tutorial 4: 01_LineBreeding/01_PhenotypicSelection/04_DoubledHaploid

**New Functions Introduced:**

**Mating:**
- `makeDH()` - Create doubled haploids (instant homozygosity)

---

## Tutorial 5: 01_LineBreeding/02_GenomicSelection

**New Functions Introduced:**

**Genomic Selection:**
- `setEBV()` - Set estimated breeding values from genomic prediction model
- `RRBLUP()` - Ridge regression BLUP model for genomic prediction

**Population Management:**
- `c()` - Combine populations (base R, used with AlphaSimR Pop objects)

**Slot Access:**
- `@nInd` - Access number of individuals slot directly

---

## Tutorial 6: 01_LineBreeding/03_TwoPartGS

**New Functions Introduced:**

*No new functions* - This tutorial uses the same functions as Tutorial 5, but demonstrates a two-part breeding strategy combining population improvement with product development.

---

## Tutorial 7: 02_ClonalBreeding/01_PhenotypicSelection

**New Functions Introduced:**

**Population Creation:**
- `runMacs2()` - Generate initial haplotypes for outcrossing/clonal species (alternative to `runMacs()`)

**Simulation Parameters:**
- `SP$addTraitADG()` - Add additive-dominance genetic trait with GxE effects (alternative to `SP$addTraitAG()`)

**Analysis:**
- `gv()` - Get genetic values as a function (alternative to accessing `@gv` slot)
- `pheno()` - Get phenotypes as a function (alternative to accessing `@pheno` slot)

**Note:** `setPheno()` is used with the `p` parameter for year/environment effects.

---

## Tutorial 8: 02_ClonalBreeding/02_PedigreeSelection

**New Functions Introduced:**

**Analysis:**
- `ebv()` - Get estimated breeding values as a function (alternative to accessing `@ebv` slot)

**Slot Access:**
- `@fixEff` - Access fixed effects slot
- `@id` - Access individual IDs slot
- `@father` - Access father IDs slot
- `@mother` - Access mother IDs slot

**Note:** This tutorial uses `setEBV()` with pedigree-based models (can use internal AlphaSimR solver or external packages like asreml).

---

## Tutorial 9: 02_ClonalBreeding/03_GenomicSelection

**New Functions Introduced:**

*No new functions* - This tutorial applies genomic selection (`RRBLUP()` and `setEBV()`) to clonal breeding, using functions already introduced in Tutorials 5 and 7-8.

---

## Tutorial 10: 03_HybridBreeding/01_PhenotypicSelection

**New Functions Introduced:**

**Simulation Parameters:**
- `SP$setVarE()` - Set permanent error variance for yield trials

**Phenotyping:**
- `setPhenoGCA()` - Set phenotypes using general combining ability (for hybrid testcross evaluation)

**Mating:**
- `hybridCross()` - Create hybrid crosses between two populations
- `runMacs()` with `split` parameter - Create heterotic pools for hybrid breeding

**Analysis:**
- `calcGCA()` - Calculate general combining ability from hybrid performance

**Slot Access:**
- `@mother` - Access mother IDs slot (used for tracking hybrid parentage)
- `@pheno` - Access phenotypes slot directly
- `@gv` - Access genetic values slot directly

**Note:** `runMacs()` is used with the `split` parameter to create separate heterotic pools, and `SP$addTraitADG()` is used (introduced in Tutorial 7) to model dominance effects important for hybrid breeding.

---

## Tutorial 11: 03_HybridBreeding/02_GenomicSelection

**New Functions Introduced:**

**Slot Access:**
- `@ebv` - Access estimated breeding values slot directly

*Note: All other functions were already introduced. This tutorial combines hybrid breeding (Tutorial 10) with genomic selection (Tutorial 5).*

---

## Tutorial 12: 03_HybridBreeding/03_TwoPartGS

**New Functions Introduced:**

*No new functions* - This tutorial demonstrates a two-part strategy for hybrid breeding with genomic selection, combining concepts from Tutorials 5, 6, and 10-11.

---

## Summary: Function Introduction Timeline

### Tutorial 1 (Foundation)
- Population creation: `runMacs()`, `newPop()`
- Simulation setup: `SimParam$new()`, `SP$restrSegSites()`, `SP$addSnpChip()`, `SP$addTraitAG()`, `SP$setTrackPed()`
- Basic operations: `randCross()`, `setPheno()`, `selectInd()`, `meanG()`, `varG()`, `cor()`

### Tutorial 2
- `self()`, `selectWithinFam()`

### Tutorial 3
- `meanP()`, `mergePops()`, `nInd()`

### Tutorial 4
- `makeDH()`

### Tutorial 5
- `setEBV()`, `RRBLUP()`, `c()` (combine populations), `@nInd`

### Tutorial 7
- `runMacs2()`, `SP$addTraitADG()`, `gv()`, `pheno()`

### Tutorial 8
- `ebv()`, `@fixEff`, `@id`, `@father`, `@mother`

### Tutorial 10
- `SP$setVarE()`, `setPhenoGCA()`, `hybridCross()`, `calcGCA()`, `@mother`, `@pheno`, `@gv`

### Tutorial 11
- `@ebv`

### Tutorials 6, 9, 12
- No new functions (demonstrate advanced strategies using previously introduced functions)

---

## Learning Path Recommendations

1. **Beginner Path:** Tutorials 1 → 2 → 3 → 4 (progressive complexity in line breeding)
2. **Genomic Selection Path:** Tutorials 1-4 → 5 (add genomic selection)
3. **Clonal Breeding Path:** Tutorials 1 → 7 → 8 → 9 (switch to clonal species)
4. **Hybrid Breeding Path:** Tutorials 1, 4 → 10 → 11 → 12 (hybrid breeding with GS)
5. **Complete Path:** All tutorials in order (comprehensive understanding)

