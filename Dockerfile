ARG BASEIMAGE=ontresearch/base-workflow-image:v0.1.1
FROM $BASEIMAGE
ARG ENVFILE=environment.yaml

COPY $ENVFILE $HOME/environment.yaml
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install -n base --file $HOME/environment.yaml \
    && micromamba clean --all --yes \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME \
    && rm -rf $CONDA_DIR/conda-meta \
    && rm -rf $CONDA_DIR/include \
    && rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip \
    && find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+'

USER $WF_UID
WORKDIR $HOME
