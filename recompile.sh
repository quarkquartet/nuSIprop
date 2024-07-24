rm -rf build/
rm nuSIprop.cpp
rm nuSIprop.cpython-310-darwin.so
python3 setup.py build_ext --inplace
