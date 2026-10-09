pragma Singleton
import QtQuick

// Thème « granite & teal » : toutes les couleurs et métriques de l'interface sont ici.
QtObject {
    // fonds
    readonly property color bg0: "#0C1114"      // fond de vue
    readonly property color bg1: "#10181C"      // fenêtre
    readonly property color bg2: "#162126"      // panneaux
    readonly property color bg3: "#1C2A30"      // survol
    readonly property color border: "#2B3B42"
    readonly property color accent: "#1FA7B5"
    readonly property color accentDoux: Qt.rgba(accent.r, accent.g, accent.b, 0.16)
    readonly property color textPrimary: "#ECF3F4"
    readonly property color textSecondary: "#9FB3B8"
    readonly property color textMuted: "#6A8086"
    readonly property color surAccent: "#06181B"          // texte posé sur l'accent
    readonly property color voile: Qt.rgba(bg0.r, bg0.g, bg0.b, 0.78)   // derrière les surimpressions
    readonly property color grille: Qt.rgba(1, 1, 1, 0.06)

    // données : modes de rupture des joints
    readonly property color modeTraction: "#E5A54B"
    readonly property color modeCisaillement: "#E77C8D"
    readonly property color modePreRompu: "#A78BFA"
    // palette arc-en-ciel style Abaqus (ancrages bleu -> rouge, découpée en bandes)
    readonly property var arcEnCiel: ["#0000FF", "#005DFF", "#00B9FF", "#00FFE8", "#00FF8B",
                                      "#00FF2E", "#2EFF00", "#8BFF00", "#E8FF00", "#FFB900",
                                      "#FF5D00", "#FF0000"]

    // typographie
    // pas de Qt.fontFamilies() ici (≈ 200 ms au démarrage) : Qt retombe seul sur Segoe UI
    readonly property string police: "Segoe UI Variable Display"
    readonly property string policeChiffres: "Segoe UI Variable Text"
    readonly property int tCorps: 13
    readonly property int tSecond: 11
    readonly property int tOverline: 10
    readonly property int tTitre: 15

    // formes et espacements (multiples de 4)
    readonly property int rCarte: 8
    readonly property int rControle: 6
    readonly property int e1: 4
    readonly property int e2: 8
    readonly property int e3: 12
    readonly property int e4: 16
    readonly property int hControle: 28
}
