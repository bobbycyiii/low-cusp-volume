FROM computop/sage:9.3
RUN git clone https://github.com/bobbycyiii/SnapPy; cd ~/SnapPy; git checkout greedy_choice
RUN cd ~/SnapPy; sage -python3 setup.py build
RUN cd ~/SnapPy; sage -python3 -m pip install .
