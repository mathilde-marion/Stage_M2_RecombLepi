# Stage M2 – LECA – Paysages de recombinaison chez les Lépidoptères : dynamiques évolutives et conséquences génomiques

## Présentation générale

Ce dépôt contient les scripts développés dans le cadre de mon stage de Master 2, portant sur l’analyse comparative de la variation et de l’évolution des taux de recombinaison chez les Lépidoptères.

Ces scripts permettent d’inférer les taux de recombinaison à l’échelle du génome, puis d’analyser leur organisation et leur évolution en lien avec différents processus génétiques et évolutifs (structure chromosomique, proximité des gènes, composition nucléotidique, efficacité de la sélection et conservation interspécifique).

L’objectif du projet est de caractériser les déterminants des paysages de recombinaison chez ce groupe, et d’évaluer dans quelle mesure ces mécanismes diffèrent de ceux décrits chez d’autres eucaryotes.

---

## Organisation du projet

Les scripts sont organisés par type d’analyse, en cohérence avec les différentes sections du matériel et méthodes du rapport.

- **bcftools_acp_smcpp**  
  Filtrage des variants (bcftools), analyses de structure des populations (ACP), inférences démographiques (SMC++) et préparation des données pour ReLERNN.

- **relernn**  
  Inférence des taux de recombinaison à l’aide de ReLERNN.

- **telomere_distance**  
  Analyse de l’organisation du taux de recombinaison le long des chromosomes (fenêtres glissantes et distance aux télomères).

- **DistanceGene_RecombRate**  
  Analyse de la variation du taux de recombinaison en fonction de la distance aux gènes, avec tests par permutation.

- **GC_content**  
  Analyse de la relation entre recombinaison et contenu en GC.

- **blast_synteny_conservation**  
  Identification de segments homologues par synténie et analyse de la conservation des taux de recombinaison entre espèces.

- **genetic_diversity_pi**  
  Estimation de la diversité génétique nucléotidique (π).

- **snpEff_piN_piS**  
  Analyse de l’efficacité de la sélection à partir du ratio πN/πS et de sa relation avec la recombinaison.

---

## Reproductibilité

Les analyses sont organisées sous forme de pipelines (Snakemake et scripts bash/R/Python), permettant une reproductibilité des résultats sur d'autres espèces, à partir des données génomiques d’entrée (VCF, GFF, FASTA).

Chaque dossier contient un README décrivant brièvement l’analyse associée et les scripts utilisés.

---

## Auteur

Mathilde Marion  
Master 2 – Biodiversité, Écologie, Évolution - Darwin
Faculté des sciences de Montpellier
