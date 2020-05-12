
def true_mag(elem: str) -> int:
    emag = dict(Ni= 2, Rb= 1, Pt= 2, Ru= 4, S= 2, Na= 1, Nb= 5, Mg= 0,
            Li= 1, Pb= 2, Pd= 0, Ti= 2, Te= 2, Rh= 3, Ta= 3, Be= 0,
            Ba= 0, As= 3, Fe= 4, Br= 1, Sr= 0, Mo= 6, He= 0, C= 2,
            B= 1, P= 3, F= 1, I= 1, H= 1, K= 1, Mn= 5, O= 2, Ne= 0,
            Kr= 0, Si= 2, Sn= 2, W= 4, V= 3, Sc= 1, N= 3, Os= 4,
            Se= 2, Zn= 0, Co= 3, Ag= 1, Cl= 1, Ca= 0, Ir= 3, Al= 1,
            Cd= 0, Ge= 2, Ar= 0, Au= 1, Zr= 2, Ga= 1, In= 1, Cs= 1,
            Cr= 6, Cu= 1, Y = 1, Sb=3,Xe=0, Hf=2, Re=5,
            Hg=0,Tl=1, La=0, Tc=5, Ce=2)
    return emag[elem]
