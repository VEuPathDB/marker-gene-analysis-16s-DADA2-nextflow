#!/usr/bin/env bash

set -euo pipefail
Rscript /usr/bin/mergeAsvsAndAssignToOtus.R \
  --asvRdsInDir . \
  --assignTaxonomyRefPath $params.trainingSet \
  --addSpeciesRefPath $params.speciesAssignment \
  --outPath ./"$genomeName"_output
