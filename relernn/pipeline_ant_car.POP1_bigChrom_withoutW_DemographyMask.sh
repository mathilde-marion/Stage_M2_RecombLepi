SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN_MOD"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="2.9e-9"
URTR="1"
DIR="./output_ant_car.POP1_bigChrom_withoutW_demographyMask/"
VCF="./ant_car.var.biallelic.fmiss80.renamed.POP1.vcf"
GENOME="./ant_car_bigChrom_withoutW.genome.bed"
MASK="./ant_car.smc_mask.bed"
DEMOGRAPHY="./Ne_ant_car.POP1.csv"

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --mask ${MASK} \
    --demographicHistory ${DEMOGRAPHY} \
    -l 1.0 \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --nTrain 100000 \
    --nVali 1000 \
    --nTest 1000 \
    --seed ${SEED} \
    -t 30

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --seed ${SEED} \
    --resume

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED}

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --nSlice 100 \
    --nReps 1000 \
    --seed ${SEED} \
    -t 30
