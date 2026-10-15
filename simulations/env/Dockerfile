# Simulation-pipeline runtime. Pure-Python numeric stack (numpy/scipy/matplotlib); no compiled extras.
# Build (from simulations/):  docker build -t smooth-simulations:latest -f env/Dockerfile .
# The repository and the FreeSurfer fsaverage5 surfaces are bind-mounted at run time, so the image
# ships only the interpreter + libraries.
FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /repo/simulations
ENV SMOOTH_REPO=/repo PROJ=/repo/simulations PYTHONPATH=/repo/simulations/code SENS_FS_SURF=/fs5surf
# run e.g. (from the repository root):
#   docker run --rm -v $PWD:/repo -v $FREESURFER_HOME/subjects/fsaverage5/surf:/fs5surf smooth-simulations \
#       python code/simulation/fwer.py --mode beta --task all
CMD ["python", "-c", "import numpy, scipy, matplotlib; print('smooth-simulations container OK', numpy.__version__)"]
