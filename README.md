# GAWM: Genome Analyst with Mongodb

## Start mongod

```bash
# rm ~/share/mongodb/data/mongod.lock
numactl --interleave=all ~/share/mongodb/bin/mongod --config ~/share/mongodb/mongod.cnf
```

## GC wave and bed count

```bash
cd ~/Scripts/gawm

# drop db if exists
mongo S288c_gc --eval "db.dropDatabase();"

# generate
perl gen_mg.pl -d S288c_gc -n S288c --dir ~/data/alignment/yeast_combine/S288C --parallel 1

# GC
perl insert_gcwave.pl -d S288c_gc --batch 1 --parallel 4

# CV
perl update_sw_cv.pl -d S288c_gc --batch 1 --parallel 4

# bed count
perl count_bed.pl -d S288c_gc --run insert -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 1 --parallel 1
perl count_bed.pl -d S288c_gc --run count --batch 1 --parallel 4

# stats
perl stat_mg.pl -d S288c_gc
```

## genome features: ofgsw

```bash
cd ~/Scripts/gawm

mongo S288c_spo11 --eval "db.dropDatabase();"

perl gen_mg.pl -d S288c_spo11 -n S288c --dir ~/data/alignment/yeast_combine/S288C  --length 10000000 --parallel 4

perl insert_bed.pl -d S288c_spo11 -tag spo11 -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 1 --parallel 4

perl stat_mg.pl -d S288c_spo11 --by tag

perl chart_mg.pl --replace ofg="DSBs" -i d:\wq\GC\autochart\131230_repli_insert\S288C_spo11.mg.xlsx
```

## Prof GC wave

```bash
cd ~/Scripts/gawm

# drop db if exists
mongo S288c_prof --eval "db.dropDatabase();"

# generate
perl gen_mg.pl -d S288c_prof -n S288c --dir ~/data/alignment/yeast_combine/S288C/chrVII.fa --parallel 1

# GC
perl -d:NYTProf insert_gcwave.pl -d S288c_prof --batch 1 --parallel 4
nytprofhtml --open

# CV
perl -d:NYTProf update_sw_cv.pl -d S288c_prof --batch 1 --parallel 4
nytprofhtml --open

# stats
perl stat_mg.pl -d S288c_prof
```
