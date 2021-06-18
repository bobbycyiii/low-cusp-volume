FROM computop/sage
RUN git clone https://github.com/bobbycyiii/SnapPy; cd ~/SnapPy; git checkout greedy_choice
RUN cd ~/SnapPy; sage -python3 setup.py build
RUN cd ~/SnapPy; sage -python3 -m pip install --user .
