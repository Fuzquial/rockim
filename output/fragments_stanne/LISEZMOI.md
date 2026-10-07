Animation des tetraedres de roche aux positions exportees (deja deformees), echelle 1.
La connectivite est reconstruite en supprimant les faces partagees dont le joint est mort (dead, sinon tBreak), car le champ fragment natif utilise D<1.
Les couleurs indiquent le signe de la vitesse verticale moyenne par tetraedre, pas une identite stable de fragment.
Un morceau sans liaison cohesive peut rester coince ou en contact : etre detache ne prouve pas un vol libre.
Le compteur au-dessus utilise tous les sommets du fragment a z>0,05 mm et le sens de sa vitesse moyenne ponderee par volume. Ce sont des indices cinematiques, pas un test de contact nul.
La surface grise est la portion restante de la surface initiale du massif, pas une reconstruction complete des faces nouvellement exposees. Insert masque pour la visibilite.
Trames figees 0-12, jusqu a 179,994 us ; aucun fichier du run modifie.
Reproduction : python tools/animate_detached_rock.py
