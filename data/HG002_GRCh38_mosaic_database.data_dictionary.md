| Field                    | Data Type | Description                                                                                      |
|--------------------------|-----------|--------------------------------------------------------------------------------------------------|
| CHROM                    | string    | Chromosome where the variant is located.                                                         |
| POS                      | integer   | Position of the variant on the chromosome.                                                       |
| REF                      | string    | Reference allele.                                                                                |
| ALT                      | string    | Alternate allele.                                                                                |
| indel_change             | string    | Change type of the variant (e.g. insertion, deletion).                                           |
| var_type                 | string    | Type of variant (e.g. snv, insertion).                                                           |
| var_cat                  | string    | Category of the variant (putative or candidate).                                                  |
| SomaticEVS               | float     | Somatic EVS value from strelka2                                                                   |
| AF                       | float     | Allele frequency of the variant.                                                                  |
| var_id                   | string    | Twinstrand Panel target IDs                                                                       |
| dup_var_hit              | string    | Variant detection with Duplex Sequencing over Twinstrand Panel target. (e.g. no_hit, hit).       |
| VAF_Category             | string    | Variant allele frequency category.                                                                 |
| SomaticEVS_Category      | string    | Somatic EVS category.                                                                             |
| Coverage_Category_Duplex | string    | Coverage category in Duplex.                                                                      |
| Genomic_Context          | string    | Genomic context of the variant.                                                                   |
| Detection_Category       | string    | Detection category of the variant.                                                                |
| easy_to_map              | string    | Easy to map region status.                                                                        |
| homopolymer_region       | string    | Homopolymer region status.                                                                        |
| germline_insertion       | string    | Germline insertion status.                                                                        |
| germline_deletion        | string    | Germline deletion status.                                                                         |
| Duplex_depth             | integer   | Number of reads at this position in Duplex.                                                       |
| Duplex_alt_depth         | integer   | Number of alternative reads at this position in Duplex.                                           |
| Duplex_mut_freq          | float     | Mutation frequency of the variant in Duplex.                                                       |
| Duplex_qual              | integer   | Quality of the variant in Duplex.                                                                 |
| Duplex_Lower_CI          | float     | Lower confidence interval in Duplex.                                                               |
| Duplex_Upper_CI          | float     | Upper confidence interval in Duplex.                                                               |
| BGI_depth                | integer   | Number of reads at this position in BGI.                                                           |
| BGI_alt_depth            | integer   | Depth of the variant the variant in BGI.                                                           |
| BGI_mut_freq             | float     | Mutation frequency of the variant in BGI.                                                           |
| BGI_map_qual             | integer   | Mapping quality of the variant in BGI.                                                             |
| BGI_base_qual            | integer   | Base quality of the variant in BGI.                                                                |
| BGI_Lower_CI             | float     | Lower confidence interval in BGI.                                                                  |
| BGI_Upper_CI             | float     | Upper confidence interval in BGI.                                                                  |
| Element_depth            | integer   | Number of reads at this position in Element.                                                       |
| Element_alt_depth        | integer   | Number of alternative reads at this position in Element.                                           |
| Element_mut_freq         | float     | Mutation frequency of the variant in Element.                                                       |
| Element_map_qual         | integer   | Mapping quality of the variant in Element.                                                         |
| Element_base_qual        | integer   | Base quality of the variant in Element.                                                            |
| Element_Lower_CI         | float     | Lower confidence interval in Element.                                                              |
| Element_Upper_CI         | float     | Upper confidence interval in Element.                                                              |
| Element_long_depth       | integer   | Number of reads at this position in long insert Element.                                           |
| Element_long_alt_depth   | integer   | Number of alternative reads at this position in long insert Element.                               |
| Element_long_mut_freq    | float     | Mutation frequency of the variant in long insert Element.                                           |
| Element_long_map_qual    | integer   | Mapping quality of the variant in long insert Element.                                              |
| Element_long_base_qual   | integer   | Base quality of the variant in long insert Element.                                                 |
| Element_long_Lower_CI    | float     | Lower confidence interval in long insert Element.                                                   |
| Element_long_Upper_CI    | float     | Upper confidence interval in long insert Element.                                                   |
| Element_standard_depth  | integer   | Number of reads at this position in standard insert Element.                                        |
| Element_standard_alt_depth | integer | Number of alternative reads at this position in standard insert Element.                            |
| Element_standard_mut_freq  | float   | Mutation frequency of the variant in standard insert Element.                                       |
| Element_standard_map_qual  | integer | Mapping quality of the variant in standard insert Element.                                          |
| Element_standard_base_qual | integer | Base quality of the variant in standard insert Element.                                             |
| Element_standard_Lower_CI  | float   | Lower confidence interval in standard insert Element.                                                |
| Element_standard_Upper_CI  | float   | Upper confidence interval in standard insert Element.                                                |
| Pacbio_depth               | integer | Number of reads at this position in Pacbio.                                                        |
| Pacbio_alt_depth           | integer | Number of alternative reads at this position in Pacbio.                                            |
| Pacbio_mut_freq            | float   | Mutation frequency of the variant in Pacbio.                                                        |
| Pacbio_map_qual            | integer | Mapping quality of the variant in Pacbio.                                                           |
| Pacbio_base_qual           | integer | Base quality of the variant in Pacbio.                                                              |
| Pacbio_Lower_CI            | float   | Lower confidence interval in Pacbio.                                                                |
| Pacbio_Upper_CI            | float   | Upper confidence interval in Pacbio.                                                                |
| Pacbio_revio_depth         | integer | Number of reads at this position in Pacbio revio.                                                  |
| Pacbio_revio_alt_depth     | integer | Number of alternative reads at this position in Pacbio revio.                                      |
| Pacbio_revio_mut_freq      | float   | Mutation frequency of the variant in Pacbio revio.                                                  |
| Pacbio_revio_map_qual      | integer | Mapping quality of the variant in Pacbio revio.                                                     |
| Pacbio_revio_base_qual     | integer | Base quality of the variant in Pacbio revio.                                                        |
| Pacbio_revio_Lower_CI      | float   | Lower confidence interval in Pacbio revio.                                                          |
| Pacbio_revio_Upper_CI      | float   | Upper confidence interval in Pacbio revio.                                                          |
| Pacbio_sequel_depth        | integer | Number of reads at this position in Pacbio sequel.                                                  |
| Pacbio_sequel_alt_depth    | integer | Number of alternative reads at this position in Pacbio sequel.                                      |
| Pacbio_sequel_mut_freq     | float   | Mutation frequency of the variant in Pacbio sequel.                                                  |
| Pacbio_sequel_map_qual     | integer | Mapping quality of the variant in Pacbio sequel.                                                     |
| Pacbio_sequel_base_qual    | integer | Base quality of the variant in Pacbio sequel.                                                        |
| Pacbio_sequel_Lower_CI     | float   | Lower confidence interval in Pacbio sequel.                                                          |
| Pacbio_sequel_Upper_CI     | float   | Upper confidence interval in Pacbio sequel.                                                          |
| Illumina_depth             | integer | Number of reads at this position in Illumina.                                                       |
| Illumina_alt_depth         | integer | Number of alternative reads at this position in Illumina.                                           |
| Illumina_mut_freq          | float   | Mutation frequency of the variant in Illumina.                                                       |
| Illumina_map_qual          | integer | Mapping quality of the variant in Illumina.                                                          |
| Illumina_base_qual         | integer | Base quality of the variant in Illumina.                                                             |
| Illumina_Lower_CI          | float   | Lower confidence interval in Illumina.                                                               |
| Illumina_Upper_CI          | float   | Upper confidence interval in Illumina.                                                               |
| Germline_Indel_Check      | string  | Germline indel check status.                                                                       |
| Duplex_Presence           | string  | Duplex presence status.                                                                            |
| Element_Presence          | string  | Element presence status.                                                                           |
| BGI_Presence              | string  | BGI presence status.                                                                               |
| Pacbio_Presence           | string  | Pacbio presence status.                                                                            |
| Illumina_Presence         | string  | Illumina presence status.                                                                          |
| x_ci                      | integer | Number of "successes" for confidence intervals.                                                     |
| n_ci                      | integer | Number of "observations" for confidence intervals.                                                   |
| x_ci_no_Dup               | integer | Number of "successes" for confidence intervals without Duplex.                                       |
| n_ci_no_Dup               | integer | Number of "observations" for confidence intervals without Duplex.                                     |
| Orthogonal_Lower_CI       | float   | Lower confidence interval from orthogonal methods.                                                   |
| Orthogonal_Upper_CI       | float   | Upper confidence interval from orthogonal methods.                                                   |
| Ortho_NO_DUP_Lower_CI     | float   | Lower confidence interval from orthogonal methods without Duplex.                                     |
| Ortho_NO_DUP_Upper_CI     | float   | Upper confidence interval from orthogonal methods without Duplex.                                     |
| Overlap_Check             | string  | Check for overlap between Orthogonal CI and Orthogonal CI without Duplex                            |
| Ortho_CI_Flag             | string  | Check whether strelka2 reported AF is within the Orthogonal CI                                      |
| NonDup_CI_Flag            | string  | Check whether strelka2 reported AF is within the Non-Duplex Orthogonal CI                           |
