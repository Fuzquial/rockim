// Joints rompus dessinés en rubans d'épaisseur constante à l'écran.
// Chaque segment = 2 triangles ; VERTEX = extrémité propre, UV0 = l'autre extrémité,
// UV1.x = côté du ruban (+1 / -1, déjà orienté pour que les deux bouts s'accordent).
VARYING vec4 vCouleur;

void MAIN()
{
    vec4 p = MODELVIEWPROJECTION_MATRIX * vec4(VERTEX, 1.0);
    vec4 q = MODELVIEWPROJECTION_MATRIX * vec4(UV0, VERTEX.z, 1.0);
    vec2 d = (q.xy / q.w - p.xy / p.w) * viewportPx;
    vec2 n = vec2(-d.y, d.x) / max(length(d), 1e-6);
    p.xy += n * UV1.x * demiLargeurPx * 2.0 / viewportPx * p.w;
    POSITION = p;
    vCouleur = COLOR;
}
