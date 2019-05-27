# nanomon-sv

## dependency

install library for htslib
```
sudo apt-get install -y zlib1g-dev libbz2-dev liblzma-dev
```

install library for samtools
```
sudo apt-get install -y libncurses5-dev`
```

 - htslib 1.9
 - samtools 1.9
 - pysam 0.15.0
 - python >= 3.5

## install 

git clone this repository

## run

0. 準備

 - minimap2 によりbamをアラインメントしておく
 - 比較する場合は genomon-sv の出力結果を別途用意しておく

1. ロングリードのbamからブレークポイントを検出する

minimap2 によりアラインメントされたbamを使用すること

```
python ./run-nanomon-sv.py parse -i ./minimap2.bam -o ./output/ 
```

2. 1. で作成されたブレークポイントと別途用意したブレークポイントを比較する

```
mkdir -p output
python ./run-nanomon-sv.py fetch --input_bp xxx.genomonSV.result.filt.txt --parsed_file xxx.junction.sort.gz --output_prefix ./output/genomon.fetch
```
