### Install prerequisites
[DGtal](https://www.dgtal.org/download/)

### Build and run

```bash
mkdir build && cd build && cmake -DGTAL_WITH_CAIRO=ON ..
make
./main
```

### For python script

```bash
python3 -m venv venv
source venv/bin/activate
pip install pandas matplotlib

python3 plot_rice.py
```

