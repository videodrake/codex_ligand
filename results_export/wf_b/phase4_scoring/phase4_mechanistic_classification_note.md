# Phase 4 Mechanistic Classification Note

## 1. Classification Logic

Each candidate site is assigned a mechanistic class based on:
- **Patch relationship** (Phase 2): orthosteric / rim / allosteric / low_relevance
- **Hotspot overlap** (Phase 1 x Phase 2): number of PPI hotspot residues in pocket
- **Docking evidence** (Phase 3): ligand support strength
- **Druggability tier** (Phase 2): for allosteric classification

## 2. Mechanistic Classes

| Class | Description | Minimum Evidence |
|-------|-------------|-----------------|
| orthosteric_disruptor_candidate | Site directly overlaps the MYO1D attachment patch. Ligand binding here would ste... | relationship_class = orthosteric_candidate AND hotspot_overl... |
| interface_rim_modulator_candidate | Site is at the rim of the PPI interface. Ligand binding may indirectly weaken MY... | relationship_class = rim_candidate AND hotspot_overlap_count... |
| allosteric_modulator_candidate | Site is spatially distant from the PPI patch but within the kinase domain. Ligan... | relationship_class = allosteric_candidate AND druggability_t... |
| ligandable_but_ppi_irrelevant_candidate | Site is druggable but has no mechanistic link to MYO1D attachment. Useful for se... | relationship_class = low_relevance_candidate AND hotspot_ove... |
| uncertain_mechanism_candidate | Evidence is insufficient or contradictory to assign a confident mechanistic labe... | Default when no other class matches... |

## 3. Classification Results

| Mechanistic Class | Count | Pockets | Confidence |
|-------------------|-------|---------|------------|
| interface_rim_modulator_candidate | 63 | 3GT8_raw_PKT01, 3GT8_raw_PKT02, 3GT8_raw_PKT05, 3GT8_raw_PKT07, 3GT8_raw_PKT08, 3GT8_raw_PKT11, 3GT8_raw_PKT14, 3GT8_raw_PKT20, EGFR_160-185_PKT02, EGFR_160-185_PKT07, EGFR_160-185_PKT09, EGFR_160-185_PKT12, EGFR_160-185_PKT15, EGFR_160-185_PKT16, EGFR_160-185_PKT23, EGFR_160-185_PKT34, EGFR_170-200_PKT06, EGFR_170-200_PKT10, EGFR_170-200_PKT12, EGFR_170-200_PKT17, EGFR_170-200_PKT28 | low, medium |
| allosteric_modulator_candidate | 6 | 3GT8_raw_PKT10, EGFR_170-200_PKT34 | low, medium |
| ligandable_but_ppi_irrelevant_candidate | 51 | 3GT8_raw_PKT04, 3GT8_raw_PKT06, 3GT8_raw_PKT09, 3GT8_raw_PKT12, 3GT8_raw_PKT13, 3GT8_raw_PKT15, 3GT8_raw_PKT17, 3GT8_raw_PKT18, 3GT8_raw_PKT19, 3GT8_raw_PKT21, EGFR_160-185_PKT06, EGFR_160-185_PKT08, EGFR_160-185_PKT10, EGFR_160-185_PKT17, EGFR_160-185_PKT18, EGFR_160-185_PKT20, EGFR_160-185_PKT21, EGFR_160-185_PKT22, EGFR_160-185_PKT24, EGFR_160-185_PKT25, EGFR_160-185_PKT26, EGFR_160-185_PKT28, EGFR_160-185_PKT31, EGFR_160-185_PKT32, EGFR_160-185_PKT33, EGFR_160-185_PKT35, EGFR_160-185_PKT36, EGFR_160-185_PKT41, EGFR_160-185_PKT42, EGFR_160-185_PKT44, EGFR_170-200_PKT02, EGFR_170-200_PKT03, EGFR_170-200_PKT07, EGFR_170-200_PKT08, EGFR_170-200_PKT09, EGFR_170-200_PKT11, EGFR_170-200_PKT14, EGFR_170-200_PKT16, EGFR_170-200_PKT18, EGFR_170-200_PKT19, EGFR_170-200_PKT20, EGFR_170-200_PKT21, EGFR_170-200_PKT22, EGFR_170-200_PKT23, EGFR_170-200_PKT24, EGFR_170-200_PKT27, EGFR_170-200_PKT29, EGFR_170-200_PKT30, EGFR_170-200_PKT32, EGFR_170-200_PKT35, EGFR_170-200_PKT37 | high |
| uncertain_mechanism_candidate | 87 | 3GT8_raw_PKT03, 3GT8_raw_PKT16, 3GT8_raw_PKT22, EGFR_160-185_PKT01, EGFR_160-185_PKT03, EGFR_160-185_PKT04, EGFR_160-185_PKT05, EGFR_160-185_PKT11, EGFR_160-185_PKT13, EGFR_160-185_PKT14, EGFR_160-185_PKT19, EGFR_160-185_PKT27, EGFR_160-185_PKT29, EGFR_160-185_PKT30, EGFR_160-185_PKT37, EGFR_160-185_PKT38, EGFR_160-185_PKT39, EGFR_160-185_PKT40, EGFR_160-185_PKT43, EGFR_170-200_PKT01, EGFR_170-200_PKT04, EGFR_170-200_PKT05, EGFR_170-200_PKT13, EGFR_170-200_PKT15, EGFR_170-200_PKT25, EGFR_170-200_PKT26, EGFR_170-200_PKT31, EGFR_170-200_PKT33, EGFR_170-200_PKT36 | low |

## 4. Per-Candidate Detail

| Pocket | Ligand | Class | Confidence | Basis |
|--------|--------|-------|------------|-------|
| 3GT8_raw_PKT01 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT02 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT03 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT05 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT07 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT08 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT10 | 173940 | allosteric_modulator_candidate | low | allosteric_relationship; tier=tier_2; ligand_support=strong |
| 3GT8_raw_PKT11 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT14 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT16 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT20 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT22 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT01 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT02 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT03 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT05 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT07 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT08 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=moderate |
| 3GT8_raw_PKT10 | 97806 | allosteric_modulator_candidate | low | allosteric_relationship; tier=tier_2; ligand_support=strong |
| 3GT8_raw_PKT11 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT14 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT16 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT20 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=moderate |
| 3GT8_raw_PKT22 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT01 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT02 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT03 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT05 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT07 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT08 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=moderate |
| 3GT8_raw_PKT10 | VAX-C12_0 | allosteric_modulator_candidate | low | allosteric_relationship; tier=tier_2; ligand_support=strong |
| 3GT8_raw_PKT11 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT14 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| 3GT8_raw_PKT16 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT20 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| 3GT8_raw_PKT22 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT01 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT02 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_160-185_PKT03 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT04 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT05 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT07 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT09 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT11 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT12 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT13 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT14 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT15 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT16 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_160-185_PKT19 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT23 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT27 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT29 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT30 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT34 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT37 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT38 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT39 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT40 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT43 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT01 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT02 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_160-185_PKT03 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT04 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT05 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT07 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT09 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=moderate |
| EGFR_160-185_PKT11 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT12 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT13 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT14 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT15 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT16 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_160-185_PKT19 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT23 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT27 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT29 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT30 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT34 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=moderate |
| EGFR_160-185_PKT37 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT38 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT39 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT40 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT43 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT01 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT02 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_160-185_PKT03 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT04 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT05 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT07 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT09 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT11 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT12 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT13 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT14 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT15 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT16 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_160-185_PKT19 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT23 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT27 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT29 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT30 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT34 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_160-185_PKT37 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT38 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT39 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT40 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_160-185_PKT43 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT01 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT04 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT05 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT06 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_170-200_PKT10 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_170-200_PKT12 | 173940 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_170-200_PKT13 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT15 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT17 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=3; ligand_support=strong |
| EGFR_170-200_PKT25 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT26 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT28 | 173940 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=moderate |
| EGFR_170-200_PKT31 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT33 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT34 | 173940 | allosteric_modulator_candidate | medium | allosteric_relationship; tier=tier_1; ligand_support=strong |
| EGFR_170-200_PKT36 | 173940 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT01 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT04 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT05 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT06 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_170-200_PKT10 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_170-200_PKT12 | 97806 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_170-200_PKT13 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT15 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT17 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=3; ligand_support=moderate |
| EGFR_170-200_PKT25 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT26 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT28 | 97806 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=weak |
| EGFR_170-200_PKT31 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT33 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT34 | 97806 | allosteric_modulator_candidate | medium | allosteric_relationship; tier=tier_1; ligand_support=strong |
| EGFR_170-200_PKT36 | 97806 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT01 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT04 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT05 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT06 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=strong |
| EGFR_170-200_PKT10 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_170-200_PKT12 | VAX-C12_0 | interface_rim_modulator_candidate | low | rim_relationship; overlap=1; ligand_support=strong |
| EGFR_170-200_PKT13 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT15 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT17 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=3; ligand_support=strong |
| EGFR_170-200_PKT25 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT26 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT28 | VAX-C12_0 | interface_rim_modulator_candidate | medium | rim_relationship; overlap=2; ligand_support=moderate |
| EGFR_170-200_PKT31 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT33 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| EGFR_170-200_PKT34 | VAX-C12_0 | allosteric_modulator_candidate | medium | allosteric_relationship; tier=tier_1; ligand_support=strong |
| EGFR_170-200_PKT36 | VAX-C12_0 | uncertain_mechanism_candidate | low | rel=allosteric_candidate; overlap=0; tier=tier_3 |
| 3GT8_raw_PKT04 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT06 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT09 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT12 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT13 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT15 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT17 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT18 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT19 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| 3GT8_raw_PKT21 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT06 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT08 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT10 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT17 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT18 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT20 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT21 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT22 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT24 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT25 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT26 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT28 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT31 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT32 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT33 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT35 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT36 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT41 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT42 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_160-185_PKT44 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT02 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT03 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT07 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT08 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT09 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT11 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT14 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT16 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT18 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT19 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT20 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT21 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT22 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT23 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT24 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT27 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT29 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT30 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT32 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT35 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |
| EGFR_170-200_PKT37 | (none) | ligandable_but_ppi_irrelevant_candidate | high | low_relevance_relationship; no_hotspot_overlap; tier=tier_3 |

## 5. Uncertainty Preservation

87 candidate(s) classified as `uncertain_mechanism_candidate`. These are retained for future evidence rather than forced into stronger labels.

---

Generated by `egfr_pipeline.phase4.mechanistic_classification`