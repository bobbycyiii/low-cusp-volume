# Summary
This repository contains enumeration and hyperbolicity code for the paper [*Hyperbolic 3-manifolds of low cusp volume*](https://arxiv.org/abs/2109.14570).

`prove_main_theorem.py` proves the main theorem, Theorem 1.1, that hyperbolic 3-manifolds with an embedded horocusp of volume at most 2.62 is obtained by filling one of the 16 manifolds listed in the following table:

------- ------- ------- ------- ------- ------- ------- ------- 
 s596    s647    s774    s776    s780    s782    s785    v2124

 v2355   v2533   v2644   v2731   v3108   v3127   v3211   v3376
------- ------- ------- ------- ------- ------- ------- ------- 
 
`remark_nonelem_nonhyp.py` proves that not every full necklace-7 manifold that embeds nonelementarily into a hyperbolic 3-manifold is itself hyperbolic. 
Search the output for `"is not"`.
This is Remark 5.24 in the final paper.

`remark_statistics.py` proves claims on counts of cusped necklace structures of various types with 4 to 7 beads in Remark 5.25 of the paper.

 `noenum_prove_gordon.py` and `prove_gordon.py` both prove Gordon's conjecture, Theorem 1.14, that the figure-8 knot exterior is the unique 1-cusped hyperbolic 3-manifold with nine or more non-hyperbolic fillings.
`noenum_prove_gordon.py` proves it assuming Theorem 1.1.
`prove_gordon.py` proves it after reproving Theorem 1.1.

`volume_bounds.py` proves the volume bounds from the paper, Theorems 1.10 and 1.12, namely that
the 14 complete non-compace hyperbolic 3-manifolds with volume at most 3.07 are

------ ------ ------ ------ ------ ------ ------ 
 m003   m004   m006   m007   m009   m010   m011

 m015   m016   m017   m019   m022   m023   m026
------ ------ ------ ------ ------ ------ ------ 

and that the closed orientable hyperbolic 3-manifolds of volume at most 1.01749 are m003(2,1), m004(5,1), and m007(3,1).

# Requirements

Clone this repository to a directory.
Then install [Docker](https://docs.docker.com/desktop/).
Then get Dunfield and Culler's SageDocker image "with the kitchen sink" by running `docker pull computop/sage`.
Then mount this repository as a local directory when running the image.
On a Linux box, the command would be

```sh
docker run -it --mount type=bind,source="$HOME/low-cusp-volume",target=/home/sage/low-cusp-volume computop/sage
```

# Running the programs

All the scripts listed above except for `volume_bounds.py` should complete under `sage`.
For instance, to prove Gordon's conjecture, run `sage prove_gordon.py` after changing directories to the local mount of this repository.

The `volume_bounds.py` script takes an argument.
To prove the first claim, run `sage volume_bounds one_cusp`.
To prove the second claim, run `sage volume_bounds closed`.
You can also prove them in one go with `sage volume_bounds one_cusp closed`.
