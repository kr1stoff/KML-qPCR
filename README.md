# KML-qPCR

## 运行

- 下载基因组

```bash
mamba -n python3.12 run poetry run \
  python -m src.kml_qpcr --genome-set-dir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes \
  download --sci-name 'Anaplasma phagocytophilum'
```

## 测试

```bash
mamba -n python3.12 run poetry run python -m tests.test_taxonomy
```
