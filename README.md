# KML-qPCR

## 运行

- 下载基因组

```bash
poetry run python -m src.kml_qpcr download \
  --sci-name 'Anaplasma phagocytophilum' \
  --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes
```

- 基因组评估

```bash
poetry python -m src.kml_qpcr assess \
  --threads 8 \
  --sci-name 'Ehrlichia chaffeensis' \
  --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes
```

## 测试

```bash
poetry run python -m tests.test_taxonomy
```
