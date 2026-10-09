import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick3D
import QtGraphs
import Proto

// Écran « Résultats » d'un triaxial 2D : barre haute, champs à gauche, courbe à droite,
// barre temporelle en bas. Les données et le rendu GPU viennent de main.py (objet `studio`).
ApplicationWindow {
    id: fenetre
    width: 1600; height: 1000
    minimumWidth: 1100; minimumHeight: 700
    visible: true
    color: Theme.bg1
    title: "rockim · Résultats — " + studio.nomRun
    font.family: Theme.police
    font.pixelSize: Theme.tCorps

    readonly property var libellesChamps: ({
        "sigmaXX": "σxx", "sigmaYY": "σyy", "sigmaXY": "τxy",
        "vonMises": "von Mises", "epsXX": "εxx", "phase": "Phase" })

    function fmt(v, d) { return Number(v).toLocaleString(Qt.locale("C"), "f", d) }

    // ---------------------------------------------------------------- composants locaux
    component Overline: Text {
        color: Theme.textMuted
        font.pixelSize: Theme.tOverline
        font.letterSpacing: 1.2
        font.weight: Font.DemiBold
        font.capitalization: Font.AllUppercase
    }

    component Carte: Rectangle {
        color: Theme.bg2
        radius: Theme.rCarte
        border.color: Theme.border
        border.width: 1
    }

    // voile semi-opaque derrière un texte posé sur la vue des champs
    component Voile: Rectangle {
        property Item cible
        x: cible.x - Theme.e2; y: cible.y - Theme.e2
        width: cible.width + 2 * Theme.e2; height: cible.height + 2 * Theme.e2
        radius: Theme.rControle
        color: Theme.voile
    }

    component Separateur: Rectangle {
        Layout.preferredWidth: 1
        Layout.preferredHeight: 24
        color: Theme.border
    }

    // pastille du sélecteur de champ (groupe segmenté)
    component Segment: Rectangle {
        id: seg
        property string cle
        property bool actif: studio.champ === cle
        implicitWidth: lib.implicitWidth + 2 * Theme.e3
        implicitHeight: Theme.hControle - 2 * Theme.e1 + 4
        radius: Theme.rControle - 2
        color: actif ? Theme.accentDoux : (zone.containsMouse ? Theme.bg3 : "transparent")
        border.width: actif ? 1 : 0
        border.color: Theme.accent
        Text {
            id: lib
            anchors.centerIn: parent
            text: fenetre.libellesChamps[seg.cle]
            color: seg.actif ? Theme.textPrimary : Theme.textSecondary
            font.pixelSize: Theme.tCorps
            font.weight: seg.actif ? Font.DemiBold : Font.Normal
        }
        MouseArea {
            id: zone
            anchors.fill: parent
            hoverEnabled: true
            cursorShape: Qt.PointingHandCursor
            onClicked: studio.champ = seg.cle
        }
    }

    // case à cocher d'un mode de joint : pastille colorée + libellé + compte
    component PuceMode: Rectangle {
        id: puce
        property int mode
        property int rang
        property string libelle
        property color teinte
        property bool coche: studio.modesVisibles[rang]
        implicitWidth: rangee.implicitWidth + 2 * Theme.e3
        implicitHeight: Theme.hControle
        radius: Theme.rControle
        color: zoneP.containsMouse ? Theme.bg3 : "transparent"
        border.color: coche ? Qt.rgba(teinte.r, teinte.g, teinte.b, 0.55) : Theme.border
        border.width: 1
        RowLayout {
            id: rangee
            anchors.centerIn: parent
            spacing: Theme.e2
            Rectangle {
                width: 12; height: 12; radius: 3
                color: puce.coche ? puce.teinte : "transparent"
                border.color: puce.teinte
                border.width: 1.5
                Text {
                    anchors.centerIn: parent
                    visible: puce.coche
                    text: "✓"
                    color: Theme.bg0
                    font.pixelSize: 9
                    font.bold: true
                }
            }
            Text {
                text: puce.libelle
                color: puce.coche ? Theme.textPrimary : Theme.textMuted
                font.pixelSize: Theme.tCorps
            }
            Text {
                text: studio.comptesModes[puce.rang]
                color: puce.coche ? puce.teinte : Theme.textMuted
                font.pixelSize: Theme.tSecond
                font.family: Theme.policeChiffres
                font.features: { "tnum": 1 }
            }
        }
        MouseArea {
            id: zoneP
            anchors.fill: parent
            hoverEnabled: true
            cursorShape: Qt.PointingHandCursor
            onClicked: studio.montrerMode(puce.mode, !puce.coche)
        }
    }

    // chiffre clé du panneau de droite
    component Chiffre: ColumnLayout {
        property string titre
        property string valeur
        property string unite
        property color teinte: Theme.textPrimary
        spacing: 2
        Overline { text: parent.titre }
        RowLayout {
            spacing: Theme.e1
            Text {
                text: parent.parent.valeur
                color: parent.parent.teinte
                font.pixelSize: 20
                font.weight: Font.DemiBold
                font.family: Theme.policeChiffres
                font.features: { "tnum": 1 }
            }
            Text {
                text: parent.parent.unite
                color: Theme.textMuted
                font.pixelSize: Theme.tSecond
                Layout.alignment: Qt.AlignBaseline
            }
        }
    }

    // ---------------------------------------------------------------- clavier
    Shortcut { sequences: ["Left"]; onActivated: studio.frame = studio.frame - 1 }
    Shortcut { sequences: ["Right"]; onActivated: studio.frame = studio.frame + 1 }
    Shortcut { sequences: ["Space"]; onActivated: lecture.running = !lecture.running }

    Timer {
        id: lecture
        interval: 120
        repeat: true
        onTriggered: studio.frame = (studio.frame + 1) % studio.nFrames
    }

    // ---------------------------------------------------------------- mise en page
    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        // ===== barre haute
        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 56
            color: Theme.bg2
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: Theme.border }

            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: Theme.e4
                anchors.rightMargin: Theme.e4
                spacing: Theme.e3

                // marque
                RowLayout {
                    spacing: Theme.e2
                    Rectangle {
                        width: 24; height: 24; radius: 6
                        color: Theme.accentDoux
                        border.color: Theme.accent
                        Rectangle { anchors.centerIn: parent; width: 8; height: 8; radius: 2; rotation: 45; color: Theme.accent }
                    }
                    ColumnLayout {
                        spacing: 0
                        Text { text: "Résultats"; color: Theme.textPrimary; font.pixelSize: Theme.tTitre; font.weight: Font.DemiBold }
                        Overline { text: "Triaxial 2D · FDEM g1" }
                    }
                }

                Separateur {}

                // sélecteur de run
                ComboBox {
                    id: choixRun
                    Layout.preferredWidth: 220
                    Layout.preferredHeight: Theme.hControle + 4
                    model: studio.runs
                    currentIndex: studio.runs.indexOf(studio.nomRun)
                    onActivated: (i) => studio.chargerRun(studio.runs[i])
                    font.pixelSize: Theme.tCorps
                    background: Rectangle {
                        radius: Theme.rControle
                        color: choixRun.hovered ? Theme.bg3 : Theme.bg1
                        border.color: choixRun.popup.visible ? Theme.accent : Theme.border
                    }
                    contentItem: RowLayout {
                        spacing: Theme.e2
                        Overline { text: "Run"; Layout.leftMargin: Theme.e3 }
                        Text {
                            Layout.fillWidth: true
                            text: choixRun.displayText
                            color: Theme.textPrimary
                            font.pixelSize: Theme.tCorps
                            elide: Text.ElideRight
                        }
                    }
                    indicator: Text {
                        x: choixRun.width - width - Theme.e3
                        anchors.verticalCenter: parent.verticalCenter
                        text: "▾"; color: Theme.textSecondary; font.pixelSize: Theme.tCorps
                    }
                    delegate: ItemDelegate {
                        width: choixRun.width
                        height: 30
                        highlighted: choixRun.highlightedIndex === index
                        contentItem: Text {
                            text: modelData
                            color: index === choixRun.currentIndex ? Theme.accent : Theme.textPrimary
                            font.pixelSize: Theme.tCorps
                            verticalAlignment: Text.AlignVCenter
                        }
                        background: Rectangle { color: parent.highlighted ? Theme.bg3 : "transparent"; radius: 4 }
                    }
                    popup: Popup {
                        y: choixRun.height + 4
                        width: choixRun.width
                        padding: Theme.e1
                        implicitHeight: contentItem.implicitHeight + 2 * Theme.e1
                        contentItem: ListView {
                            implicitHeight: contentHeight
                            model: choixRun.delegateModel
                            clip: true
                        }
                        background: Rectangle { color: Theme.bg2; radius: Theme.rControle; border.color: Theme.border }
                    }
                }

                Separateur {}

                // sélecteur de champ (groupe segmenté)
                Overline { text: "Champ" }
                Rectangle {
                    implicitWidth: segs.implicitWidth + 2 * Theme.e1
                    implicitHeight: Theme.hControle + 4
                    radius: Theme.rControle
                    color: Theme.bg1
                    border.color: Theme.border
                    RowLayout {
                        id: segs
                        anchors.centerIn: parent
                        spacing: 2
                        Repeater {
                            model: ["sigmaXX", "sigmaYY", "sigmaXY", "vonMises", "epsXX", "phase"]
                            Segment { cle: modelData }
                        }
                    }
                }

                Separateur {}

                Overline { text: "Joints rompus" }
                PuceMode { mode: 1; rang: 0; libelle: "Traction"; teinte: Theme.modeTraction }
                PuceMode { mode: 2; rang: 1; libelle: "Cisaillement"; teinte: Theme.modeCisaillement }
                PuceMode { mode: 4; rang: 2; libelle: "Pré-rompu"; teinte: Theme.modePreRompu }

                Item { Layout.fillWidth: true }

                // badges d'identité du run
                Rectangle {
                    implicitWidth: badge.implicitWidth + 2 * Theme.e3
                    implicitHeight: 24
                    radius: 12
                    color: Theme.accentDoux
                    Text {
                        id: badge
                        anchors.centerIn: parent
                        text: "σ₃ = " + fenetre.fmt(studio.sigma3, 0) + " MPa"
                        color: Theme.accent
                        font.pixelSize: Theme.tSecond
                        font.weight: Font.DemiBold
                    }
                }
            }
        }

        // ===== zone centrale
        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.margins: Theme.e3
            spacing: Theme.e3

            // ---- vue des champs
            Carte {
                Layout.fillWidth: true
                Layout.fillHeight: true
                color: Theme.bg0
                clip: true

                Item {
                    id: vue
                    objectName: "vueChamps"
                    anchors.fill: parent
                    anchors.margins: 1

                    // état caméra : zoom relatif au cadrage, décalage en mm
                    property real zoom: 1.0
                    property real decalX: 0
                    property real decalY: 0
                    readonly property var b: studio.boite
                    readonly property real cx: (b[0] + b[2]) / 2
                    readonly property real cy: (b[1] + b[3]) / 2
                    // cadrage : la boîte occupe 88 % de la hauteur utile (barre de couleur à droite)
                    readonly property real magCadre: Math.min((width - 220) / (b[2] - b[0]),
                                                              0.88 * height / (b[3] - b[1]))
                    readonly property real mag: Math.max(1e-3, magCadre * zoom)

                    function recadrer() { zoom = 1; decalX = 0; decalY = 0 }

                    View3D {
                        id: v3d
                        anchors.fill: parent
                        renderMode: View3D.Offscreen
                        camera: camera
                        environment: SceneEnvironment {
                            backgroundMode: SceneEnvironment.Transparent
                            antialiasingMode: SceneEnvironment.MSAA
                            antialiasingQuality: SceneEnvironment.High
                            tonemapMode: SceneEnvironment.TonemapModeNone
                        }
                        OrthographicCamera {
                            id: camera
                            // le centre de la boîte est décalé vers la gauche pour laisser la barre de couleur
                            x: vue.cx + vue.decalX + 50 / vue.mag
                            y: vue.cy + vue.decalY
                            z: 500
                            clipNear: 1
                            clipFar: 1000
                            horizontalMagnification: vue.mag
                            verticalMagnification: vue.mag
                        }
                        Model {
                            geometry: GeometrieChamps { id: geoChamps }
                            materials: DefaultMaterial {
                                lighting: DefaultMaterial.NoLighting
                                vertexColorsEnabled: true
                                cullMode: Material.NoCulling
                            }
                        }
                        Model {
                            geometry: GeometrieJoints { id: geoJoints }
                            materials: CustomMaterial {
                                shadingMode: CustomMaterial.Unshaded
                                cullMode: Material.NoCulling
                                vertexShader: "joints.vert"
                                fragmentShader: "joints.frag"
                                property size viewportPx: Qt.size(v3d.width, v3d.height)
                                property real demiLargeurPx: 1.1
                            }
                        }
                    }

                    // interaction : molette = zoom centré, glisser = déplacer, double-clic = recadrer
                    MouseArea {
                        anchors.fill: parent
                        acceptedButtons: Qt.LeftButton | Qt.MiddleButton
                        cursorShape: pressed ? Qt.ClosedHandCursor : Qt.OpenHandCursor
                        property real px; property real py
                        onPressed: (m) => { px = m.x; py = m.y }
                        onPositionChanged: (m) => {
                            vue.decalX -= (m.x - px) / vue.mag
                            vue.decalY += (m.y - py) / vue.mag
                            px = m.x; py = m.y
                        }
                        onDoubleClicked: vue.recadrer()
                        onWheel: (w) => {
                            const f = Math.pow(1.0015, w.angleDelta.y)
                            const z = Math.max(0.5, Math.min(60, vue.zoom * f))
                            const m0 = vue.mag, m1 = vue.magCadre * z
                            // point monde sous la souris, conservé après le zoom
                            const dx = w.x - vue.width / 2, dy = w.y - vue.height / 2
                            vue.decalX += dx / m0 - dx / m1 + 50 / m0 - 50 / m1
                            vue.decalY -= dy / m0 - dy / m1
                            vue.zoom = z
                        }
                    }
                }

                // titre de la vue
                Voile { cible: titreVue }
                ColumnLayout {
                    id: titreVue
                    anchors { left: parent.left; top: parent.top; margins: Theme.e4 }
                    spacing: 2
                    Overline { text: "Champ par élément" }
                    Text {
                        text: fenetre.libellesChamps[studio.champ]
                              + (studio.echelleCouleurs[2] ? "  (" + studio.echelleCouleurs[2] + ")" : "")
                        color: Theme.textPrimary
                        font.pixelSize: 18
                        font.weight: Font.DemiBold
                    }
                    Text {
                        text: studio.nTri.toLocaleString(Qt.locale("fr_FR"), "f", 0) + " triangles · "
                              + studio.nJoint.toLocaleString(Qt.locale("fr_FR"), "f", 0) + " joints"
                        color: Theme.textMuted
                        font.pixelSize: Theme.tSecond
                    }
                    Text {
                        text: "t = " + fenetre.fmt(studio.tempsMs, 2) + " ms · frame " + studio.frame
                        color: Theme.textSecondary
                        font.pixelSize: Theme.tSecond
                        font.family: Theme.policeChiffres
                    }
                }

                // aide
                Voile { cible: aide }
                Text {
                    id: aide
                    anchors { left: parent.left; bottom: parent.bottom; margins: Theme.e4 }
                    text: "Molette : zoom · glisser : déplacer · double-clic : recadrer · ← → : frames"
                    color: Theme.textMuted
                    font.pixelSize: Theme.tSecond
                }

                // échelle graphique (10 mm)
                Column {
                    anchors { right: barre.left; bottom: parent.bottom; margins: Theme.e4 }
                    anchors.rightMargin: Theme.e4 * 2
                    spacing: Theme.e1
                    Text {
                        anchors.horizontalCenter: parent.horizontalCenter
                        text: "10 mm"; color: Theme.textSecondary; font.pixelSize: Theme.tSecond
                    }
                    Rectangle {
                        width: 10 * vue.mag; height: 3; radius: 1.5
                        color: Theme.textSecondary
                    }
                }

                // barre de couleur graduée
                Rectangle {
                    x: barre.x - Theme.e3; y: barre.y - Theme.e3 - 24
                    width: barre.width + Theme.e3; height: barre.height + 2 * Theme.e3 + 24
                    radius: Theme.rCarte
                    color: Theme.voile
                }
                Item {
                    id: barre
                    anchors { right: parent.right; verticalCenter: parent.verticalCenter; rightMargin: Theme.e4 }
                    width: 96
                    height: Math.min(parent.height * 0.62, 460)
                    readonly property var ech: studio.echelleCouleurs
                    readonly property int n: studio.nbBandes
                    readonly property int dec: {
                        const pas = Math.abs(ech[1] - ech[0]) / n
                        return pas >= 10 ? 0 : pas >= 1 ? 1 : pas >= 0.1 ? 2 : 3
                    }
                    Column {
                        id: bandes
                        width: 16
                        height: parent.height
                        Repeater {
                            model: barre.n
                            Rectangle {
                                width: 16
                                height: bandes.height / barre.n
                                color: Theme.arcEnCiel[barre.n - 1 - index]
                            }
                        }
                    }
                    Rectangle {
                        width: 16; height: parent.height
                        color: "transparent"; border.color: Theme.border; radius: 2
                    }
                    Repeater {
                        model: barre.n + 1
                        Item {
                            y: index * barre.height / barre.n
                            Rectangle { x: 16; width: 4; height: 1; color: Theme.textMuted }
                            Text {
                                x: 24; y: -height / 2
                                text: studio.champ === "phase"
                                      ? "" : fenetre.fmt(barre.ech[1] - index * (barre.ech[1] - barre.ech[0]) / barre.n, barre.dec)
                                color: Theme.textSecondary
                                font.pixelSize: Theme.tSecond
                                font.family: Theme.policeChiffres
                                font.features: { "tnum": 1 }
                            }
                        }
                    }
                    Overline {
                        anchors.bottom: parent.top
                        anchors.bottomMargin: Theme.e2
                        text: studio.echelleCouleurs[2] || "phase"
                    }
                }
            }

            // ---- panneau de droite : courbe + chiffres
            ColumnLayout {
                Layout.preferredWidth: 500
                Layout.maximumWidth: 500
                Layout.fillWidth: false
                Layout.fillHeight: true
                spacing: Theme.e3

                Carte {
                    Layout.fillWidth: true
                    Layout.fillHeight: true

                    ColumnLayout {
                        anchors.fill: parent
                        anchors.margins: Theme.e4
                        spacing: Theme.e2

                        RowLayout {
                            Layout.fillWidth: true
                            ColumnLayout {
                                spacing: 2
                                Overline { text: "Courbe de chargement" }
                                Text { text: "q en fonction de ε axial"; color: Theme.textPrimary; font.pixelSize: Theme.tTitre; font.weight: Font.DemiBold }
                            }
                            Item { Layout.fillWidth: true }
                            Text {
                                text: "clic sur la courbe : aller à la frame"
                                color: Theme.textMuted
                                font.pixelSize: Theme.tSecond
                            }
                        }

                        Item {
                          Layout.fillWidth: true
                          Layout.fillHeight: true
                          GraphsView {
                            id: graphe
                            anchors.fill: parent
                            marginLeft: 4; marginRight: 12; marginTop: 8; marginBottom: 4
                            readonly property var cb: studio.courbeBornes
                            theme: GraphsTheme {
                                colorScheme: GraphsTheme.ColorScheme.Dark
                                backgroundVisible: false
                                plotAreaBackgroundVisible: false
                                gridVisible: true
                                grid.mainColor: Theme.grille
                                grid.subColor: "transparent"
                                grid.mainWidth: 1
                                axisX.mainColor: Theme.border
                                axisY.mainColor: Theme.border
                                labelTextColor: Theme.textSecondary
                                labelBackgroundVisible: false
                                axisXLabelFont.pixelSize: Theme.tSecond
                                axisYLabelFont.pixelSize: Theme.tSecond
                                axisXLabelFont.family: Theme.policeChiffres
                                axisYLabelFont.family: Theme.policeChiffres
                                seriesColors: [Theme.accent]
                            }
                            axisX: ValueAxis {
                                min: 0
                                max: graphe.cb[0] * 1.04
                                tickInterval: max > 0.8 ? 0.2 : 0.1
                                labelDecimals: 2
                                subTickCount: 0
                                titleText: "ε axial (%)"
                                titleColor: Theme.textSecondary
                                titleFont.pixelSize: Theme.tSecond
                                subGridVisible: false
                            }
                            axisY: ValueAxis {
                                min: graphe.cb[1]
                                max: Math.ceil(graphe.cb[2] * 1.1 / 20) * 20
                                tickInterval: 20
                                labelDecimals: 0
                                subTickCount: 0
                                titleText: "q (MPa)"
                                titleColor: Theme.textSecondary
                                titleFont.pixelSize: Theme.tSecond
                                subGridVisible: false
                            }
                            LineSeries {
                                id: serie
                                width: 2
                                color: Theme.accent
                            }
                          }

                            // marqueur de la frame courante (repère du graphe -> pixels)
                            Item {
                                id: repere
                                readonly property rect pa: graphe.plotArea
                                readonly property real mx: pa.x + (studio.epsCourant - graphe.axisX.min)
                                                          / (graphe.axisX.max - graphe.axisX.min) * pa.width
                                readonly property real my: pa.y + pa.height - (studio.qCourant - graphe.axisY.min)
                                                          / (graphe.axisY.max - graphe.axisY.min) * pa.height
                                anchors.fill: parent
                                Rectangle {
                                    x: repere.mx; y: repere.pa.y
                                    width: 1; height: repere.pa.height
                                    color: Theme.accent; opacity: 0.35
                                }
                                Rectangle {
                                    x: repere.mx - 7; y: repere.my - 7
                                    width: 14; height: 14; radius: 7
                                    color: Theme.bg2
                                    border.color: Theme.accent; border.width: 2
                                    Rectangle { anchors.centerIn: parent; width: 6; height: 6; radius: 3; color: Theme.textPrimary }
                                }
                            }
                            MouseArea {
                                anchors.fill: parent
                                cursorShape: Qt.CrossCursor
                                function aller(x) {
                                    const pa = graphe.plotArea
                                    const e = graphe.axisX.min + (x - pa.x) / pa.width * (graphe.axisX.max - graphe.axisX.min)
                                    studio.frame = studio.frameProche(e)
                                }
                                onPressed: (m) => aller(m.x)
                                onPositionChanged: (m) => { if (pressed) aller(m.x) }
                            }
                        }
                    }
                }

                // chiffres clés
                Carte {
                    Layout.fillWidth: true
                    Layout.preferredHeight: 156
                    GridLayout {
                        anchors.fill: parent
                        anchors.margins: Theme.e4
                        columns: 3
                        rowSpacing: Theme.e3
                        columnSpacing: Theme.e4
                        Chiffre { titre: "q à la frame"; valeur: fenetre.fmt(studio.qCourant, 1); unite: "MPa"; teinte: Theme.accent }
                        Chiffre { titre: "Déformation axiale"; valeur: fenetre.fmt(studio.epsCourant, 3); unite: "%" }
                        Chiffre { titre: "q pic"; valeur: fenetre.fmt(studio.qPic, 1); unite: "MPa" }
                        Chiffre { titre: "Traction"; valeur: studio.comptesModes[0]; unite: "joints"; teinte: Theme.modeTraction }
                        Chiffre { titre: "Cisaillement"; valeur: studio.comptesModes[1]; unite: "joints"; teinte: Theme.modeCisaillement }
                        Chiffre { titre: "Pré-rompus"; valeur: studio.comptesModes[2]; unite: "joints"; teinte: Theme.modePreRompu }
                    }
                }
            }
        }

        // ===== barre temporelle
        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 64
            color: Theme.bg2
            Rectangle { anchors.top: parent.top; width: parent.width; height: 1; color: Theme.border }

            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: Theme.e4
                anchors.rightMargin: Theme.e4
                spacing: Theme.e4

                // lecture / pause
                Rectangle {
                    width: 36; height: 36; radius: 18
                    color: zoneL.containsMouse ? Qt.lighter(Theme.accent, 1.15) : Theme.accent
                    Text {
                        anchors.centerIn: parent
                        anchors.horizontalCenterOffset: lecture.running ? 0 : 1
                        text: lecture.running ? "❚❚" : "▶"
                        color: Theme.surAccent
                        font.pixelSize: lecture.running ? 11 : 14
                    }
                    MouseArea {
                        id: zoneL
                        anchors.fill: parent
                        hoverEnabled: true
                        cursorShape: Qt.PointingHandCursor
                        onClicked: lecture.running = !lecture.running
                    }
                }

                ColumnLayout {
                    spacing: 0
                    Layout.preferredWidth: 120
                    Overline { text: "Temps" }
                    Text {
                        text: fenetre.fmt(studio.tempsMs, 2) + " ms"
                        color: Theme.textPrimary
                        font.pixelSize: 18
                        font.weight: Font.DemiBold
                        font.family: Theme.policeChiffres
                        font.features: { "tnum": 1 }
                    }
                }

                // curseur temporel avec une graduation par frame
                Slider {
                    id: curseur
                    Layout.fillWidth: true
                    from: 0; to: studio.nFrames - 1; stepSize: 1
                    snapMode: Slider.SnapAlways
                    value: studio.frame
                    onMoved: studio.frame = Math.round(value)
                    focusPolicy: Qt.NoFocus
                    background: Item {
                        x: curseur.leftPadding
                        y: curseur.topPadding + curseur.availableHeight / 2 - 2
                        width: curseur.availableWidth
                        height: 4
                        Rectangle { anchors.fill: parent; radius: 2; color: Theme.bg3 }
                        Rectangle {
                            width: curseur.visualPosition * parent.width
                            height: parent.height; radius: 2
                            color: Theme.accent
                        }
                        Repeater {
                            model: studio.nFrames
                            Rectangle {
                                x: index / (studio.nFrames - 1) * parent.width - 0.5
                                y: 10
                                width: 1; height: index % 5 === 0 ? 6 : 3
                                color: index <= studio.frame ? Theme.accent : Theme.textMuted
                                opacity: 0.8
                            }
                        }
                    }
                    handle: Rectangle {
                        x: curseur.leftPadding + curseur.visualPosition * (curseur.availableWidth - width)
                        y: curseur.topPadding + curseur.availableHeight / 2 - height / 2
                        width: 16; height: 16; radius: 8
                        color: Theme.textPrimary
                        border.color: Theme.accent
                        border.width: 3
                    }
                }

                ColumnLayout {
                    spacing: 0
                    Layout.preferredWidth: 110
                    Overline { text: "Frame"; Layout.alignment: Qt.AlignRight }
                    Text {
                        Layout.alignment: Qt.AlignRight
                        text: studio.frame + " / " + (studio.nFrames - 1)
                        color: Theme.textPrimary
                        font.pixelSize: 18
                        font.weight: Font.DemiBold
                        font.family: Theme.policeChiffres
                        font.features: { "tnum": 1 }
                    }
                }
            }
        }
    }

    Component.onCompleted: studio.demarrer(geoChamps, geoJoints, serie, Theme.arcEnCiel,
                                           Theme.modeTraction, Theme.modeCisaillement, Theme.modePreRompu)
}
