# FermiHarmonics

FermiHarmonics is a Julia code for solving a toy linearized 2D Boltzmann equation for Fermi-liquid transport, including momentum-relaxing and momentum-conserving collisions, on unstructured meshes.  It is mostly a lightweight wrapper for `Trixi.jl`, which is a library for high-order PDE solutions using Discontinous Galerkin Spectral Element Method (DGSEM).  

![FermiHarmonics regime demo](demo/streamlines.png)

Documentation: <https://docs.jackhfarrell.com>

## Citation
If you happen to find this code useful, it would be great if you would cite our upcoming theory/numerics paper as well as the codebase itself!

- Code release (Zenodo): `TBD`
- Theory/application paper: `TBD`

### BibTeX (Code)

```bibtex
@software{fermiharmonics_zenodo,
  author = {Farrell, Jack H.},
  title = {FermiHarmonics},
  year = {TBD},
  publisher = {Zenodo},
  doi = {TBD},
  url = {TBD}
}
```

### BibTeX (Theory Paper)

```bibtex
@article{fermiharmonics_theory,
  author = {Farrell, Jack H. and ...},
  title = {TBD},
  journal = {TBD},
  year = {TBD},
  doi = {TBD},
  url = {TBD}
}
```
