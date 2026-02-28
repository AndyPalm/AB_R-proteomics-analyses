# Mouse-Human Protein Mapping

## Overview
This analysis maps mouse protein identifiers to their human orthologs using the JAX (Jackson Laboratory) mouse-human ortholog database. It handles cases where direct matches are not found and provides diagnostic information on mapping success.

## Inputs
- **Mouse gene list** (CSV): Mouse protein IDs and gene names from proteomics data
- **Ortholog lookup table** (XLSX): Mouse-human ortholog mapping from JAX database
- **JAX database** (online): Mouse-human sequence homology data (downloaded from informatics.jax.org)

## Methods
1. **Database Download**: Retrieves current mouse-human ortholog table from JAX website
2. **Protein Matching**: Searches for mouse protein IDs in lookup table
3. **Ortholog Retrieval**: Extracts human ortholog information for matched proteins
4. **Fallback Handling**: Preserves original mouse gene names if human ortholog not found
5. **Diagnostic Analysis**: Identifies multiple matches and mapping failures

## Key Parameters
- **Lookup key**: ProteinID (mouse UniProt accession)
- **Ortholog class**: DB.Class.Key for grouping orthologs
- **Match strategy**: First match selected if multiple orthologs exist
- **Fallback**: Original mouse gene name used if no human match found

## Outputs
- Mapped dataset with human ProteinID and gene names
- Mapping success/failure flags
- Diagnostic report of multiple matches
- CSV export of mouse-to-human mapping
