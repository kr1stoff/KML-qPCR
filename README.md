# KML-qPCR

## 总体流程

1. 下载基因组
2. 注释基因组
3. 评估基因组
4. 保守基因获取

```bash
# 1. 下载基因组
poetry run python -m src.kml_qpcr download \
  --sci-name 'Coxiella Burnetii' \
  --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes

# [可选] 1. 导入客户基因组
poetry run python -m src.kml_qpcr load \
  --threads 32 \
  --sci-name 'Chlamydia psittaci' \
  --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes \
  --customer-genome-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Chlamydia_psittaci/customer

# 2. 注释基因组
poetry run python -m src.kml_qpcr annotate \
  --sci-name 'Coxiella Burnetii' \
  --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes \
  --threads 32

# 3. 评估基因组
poetry run python -m src.kml_qpcr assess \
  --threads 8 \
  --sci-name 'Coxiella Burnetii' \
  --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes
```

## 测试

```bash
poetry run python -m tests.test_taxonomy
```
