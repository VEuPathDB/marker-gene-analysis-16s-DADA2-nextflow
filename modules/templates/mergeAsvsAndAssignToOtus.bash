#!/usr/bin/env bash

Rscript /usr/bin/mergeAsvsAndAssignToOtus.R \
  --asvRdsInDir . \
  --assignTaxonomyRefPath $params.trainingSet \
  --addSpeciesRefPath $params.speciesAssignment \
  --outPath ./"$genomeName"_output
