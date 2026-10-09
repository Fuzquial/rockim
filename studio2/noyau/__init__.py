"""Noyau de l'interface graphique v2 (spec 007) : toute la logique, aucune interface.

essai          modèle d'un essai triaxial 2D -> deck .cfg
materiaux      préréglages Red Bohus (une ou trois phases)
geometrie      plans de discontinuités, éléments par grain
validation     garde-fous du solveur et règles maison, avant lancement
maillage       maillage Gmsh non structuré (régime homogène)
file           file de calculs persistante, exécution parallèle
resultats      lecture d'un run, cache binaire pour l'affichage
depouillement  grandeurs de l'essai (q pic, E sécant, modes de rupture...)
"""
