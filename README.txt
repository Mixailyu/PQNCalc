The result of this example is the vibrational energy of the transition  [000->100] for ^18 O_3 molecule. The molecular parameters were taken from [A. Barbe, C. Secroun, and P. Jouve. Infrared spectra of ^16o3 and ^18o3: Darling anddennison resonance and anharmonic potential function of ozone.Journal of MolecularSpectroscopy, 49(2):171–182, 1974].

It is worth mentioning that the relations between different types of constants. Constants used in our calculations have the relation with analogous constants from [A. Barbe, C. Secroun, and P. Jouve. Infrared spectra of ^16o3 and ^18o3: Darling anddennison resonance and anharmonic potential function of ozone.Journal of MolecularSpectroscopy, 49(2):171–182, 1974], [A. Barbe, A. Chichery, T. Cours, Vl. G. Tyuterev, and J. J. Plateaux. Update ofthe anharmonic force field parameters of the ozone molecule.Journal of MolecularStructure, 616(1-3):55–65, 2002], [P. Hennig and G. Strey.  Anharmonic force field and isotopic relations of ozone.Zeitschrift f ̈ur Naturforschung A, 31(3-4):244–250, 1976]:


a_{j_1,j_2, ..., j_r}=k_{j_1,j_2, ..., j_r} * 2^(-(p+2)/2)=Phi_{j_1,j_2, ..., j_r} * 2^(-(p+2)/2)/(j_1! j_2! ... j_r!)


After starting the script, we suggest to follow the next scheme:

1. E0=AE([0,0,0],[0,0,0],0,'s')

2. E=AE([0,0,0],[0,0,0],2,'s')

3. E0=sy.sympify(E0)

4. E=sy.sympify(E)

5. zamena=[(n__i, 0), (n__j, 0), (n__k, 0),(omega__i, 1069.9), (omega__j, 674.9), (omega__k, 1026.8), (A__iii, -15.556), (A__ijj, -8.202), (A__ikk, -73.079), (A__iij, -9.617), (A__jjj, -6.222), (A__jkk, -19.162), (A__iiii, 0.5), (A__iijj, -0.4), (A__jjjj, 0.125), (A__iikk, 6.3), (A__jjkk, -1.3), (A__kkkk, 1.5), (A__iik,0), (A__jjk,0), (A__kkk,0), (A__ijk,0), (A__iiij, 0), (A__ijjj, 0), (A__ijkk, 0), (A__iiik, 0), (A__ijjk, 0), (A__ikkk, 0)]

6. zamena1=[(n__i, 1), (n__j, 0), (n__k, 0),(omega__i, 1069.9), (omega__j, 674.9), (omega__k, 1026.8), (A__iii, -15.556), (A__ijj, -8.202), (A__ikk, -73.079), (A__iij, -9.617), (A__jjj, -6.222), (A__jkk, -19.162), (A__iiii, 0.5),  (A__iijj, -0.4), (A__jjjj, 0.125), (A__iikk, 6.3), (A__jjkk, -1.3),  (A__kkkk, 1.5), (A__iik,0), (A__jjk,0), (A__kkk,0), (A__ijk,0), (A__iiij, 0), (A__ijjj, 0), (A__ijkk, 0), (A__iiik, 0), (A__ijjk, 0), (A__ikkk, 0)]

7. E0.subs(zamena1)+E.subs(zamena1)-(E0.subs(zamena)+E.subs(zamena))

8. OUT: 1041.74984600770


In order to calculate resonance energy levels, we suggest next set of commands for calculation of [200] and [002] resonance levels for $^{18}O_3$ molecule:


1. zamena=[(n__i, 0), (n__j, 0), (n__k, 0),(omega__i, 1069.9), (omega__j, 674.9), (omega__k, 1026.8), (A__iii, -15.556), (A__ijj, -8.202), (A__ikk, -73.079), (A__iij, -9.617), (A__jjj, -6.222), (A__jkk, -19.162), (A__iiii, 0.5), (A__iijj, -0.4), (A__jjjj, 0.125), (A__iikk, 6.3), (A__jjkk, -1.3), (A__kkkk, 1.5), (A__iik,0), (A__jjk,0), (A__kkk,0), (A__ijk,0), (A__iiij, 0), (A__ijjj, 0), (A__ijkk, 0), (A__iiik, 0), (A__ijjk, 0), (A__ikkk, 0)]

2. Energy_Resonance([2,0,0],[0,0,2], zamena)

3. OUT: (2079.553, 1946.618)


Description for each line:

1. zamena is a list containing molecular parameters.

2. Calculation of resonance energy levels by function Energy_Resonance.

3. The result of calculations as a list with two elements.