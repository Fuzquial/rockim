#pragma once
// =====================================================================
//  JointTsl.hpp — noyau de la loi cohesive EXTRINSEQUE de la note
//  « Lois constitutives proposees pour un FDEM hybride a insertion
//    adaptative », note de travail, septembre 2026, sections 2.2 et 2.4.
//
//  Fonctions PURES, sans etat, sans dependance a la dimension : le
//  solveur 2D et le solveur 3D appellent EXACTEMENT le meme code, de
//  sorte que la loi ne peut pas diverger entre les deux par recopie.
//
//  Rien ici n'est branche par defaut. Le cablage se fait derriere les
//  cles opt-in `insertionCriterion = elliptic`, `jointTSL = camacho` et
//  `jointMixLaw = bk` (principe VIII : croissance par addition).
//
//  CONVENTION DE SIGNE — traction NORMALE POSITIVE a l'ouverture,
//  conformement a la section 2 de la note (et non a la convention
//  geomecanique de sa section 1).
//
//  Unites SI strictes, comme tout rockim : m, Pa, J/m^2, s.
// =====================================================================
#include <algorithm>
#include <cmath>

namespace rockim {
namespace jtsl {

// ---------------------------------------------------------------------
//  Crochet de Macaulay <x> = max(x, 0)
// ---------------------------------------------------------------------
inline double mac(double x) { return x > 0.0 ? x : 0.0; }

// ---------------------------------------------------------------------
//  §2.2 eq. 12 — CRITERE D'INSERTION ELLIPTIQUE
//
//      Phi_F = ( <t_n> / t_n0 )^2 + ( t_s / t_s0 )^2   >= 1
//
//  C'est une NORME EFFECTIVE, pas un OU logique : une facette chargee a
//  80 % de sa resistance en traction ET a 80 % en cisaillement s'insere
//  (Phi = 1,28), la ou le OU la laisserait intacte. C'est la difference
//  de selectivite relevee dans biblio_insertion/ entre les fondateurs
//  (Camacho-Ortiz 1996) et la litterature recente.
//
//  Retourne Phi_F. Un seuil nul ou negatif est traite comme « cette
//  branche ne peut pas rompre » (contribution nulle) plutot que comme
//  une division par zero : c'est le comportement attendu d'un joint a
//  resistance infinie (bulkFt = 1e12 et consorts).
// ---------------------------------------------------------------------
inline double phiInsert(double tn, double ts, double tn0, double ts0)
{
    double rn = (tn0 > 0.0) ? mac(tn) / tn0 : 0.0;
    double rs = (ts0 > 0.0) ? std::fabs(ts) / ts0 : 0.0;
    return rn * rn + rs * rs;
}

// ---------------------------------------------------------------------
//  §2.2 eq. 13 — FACTEUR DYNAMIQUE
//
//      DIF_k(x) = max[ 1, (x / edot0)^a_k ]
//
//  avec a_s < a_t (le cisaillement est moins sensible a la vitesse que
//  la traction). Evalue UNE SEULE FOIS, a l'insertion, puis gele : le
//  gel est la responsabilite de l'appelant (champ difStamped du joint),
//  pas de cette fonction.
// ---------------------------------------------------------------------
inline double dif(double rate, double rate0, double expo)
{
    if (!(rate > 0.0) || !(rate0 > 0.0) || !(expo > 0.0)) return 1.0;
    double d = std::pow(rate / rate0, expo);
    return d > 1.0 ? d : 1.0;
}

// ---------------------------------------------------------------------
//  §2.4 eq. 16 — SEPARATION ET TRACTION EFFECTIVES (Camacho-Ortiz)
//
//      delta_m = sqrt( <delta_n>^2 + beta^2 delta_s^2 )
//      t_m     = sqrt( <t_n>^2     + beta^-2 t_s^2   )
//
//  beta = (t_s0 / t_n0) evalue A L'INSERTION (eq. 16, indice
//  « insertion »), donc fige : c'est un rapport de RESISTANCES, pas de
//  tractions courantes.
// ---------------------------------------------------------------------
inline double effSep(double dn, double ds, double beta)
{
    double a = mac(dn);
    return std::sqrt(a * a + beta * beta * ds * ds);
}

inline double effTrac(double tn, double ts, double beta)
{
    double a = mac(tn);
    double b = (beta > 0.0) ? ts / beta : 0.0;
    return std::sqrt(a * a + b * b);
}

// ---------------------------------------------------------------------
//  §2.4 eq. 18 — MIXITE DE MODE, MESUREE SUR LES TRACTIONS A L'INSERTION
//
//  L'eq. 18 definit la mixite sur les SEPARATIONS :
//      m = beta^2 delta_s^2 / ( <delta_n>^2 + beta^2 delta_s^2 )
//  mais a l'instant de l'insertion le joint nait a separation NULLE
//  (§2.3 : « continuite des vitesses et positions … aucun saut
//  cinematique »), donc m y est 0/0.
//
//  On l'evalue sur les tractions, ce qui est la MEME grandeur prise a
//  la limite de l'ouverture naissante. En effet l'eq. 19 donne
//      t_n = (t_m/delta_m) <delta_n>   et   t_s = (t_m/delta_m) beta^2 delta_s
//  donc  r = t_s/t_n = beta^2 delta_s / <delta_n>, d'ou beta^2 delta_s^2
//  = r^2 <delta_n>^2 / beta^2 et, en reportant dans l'eq. 18 :
//
//      m = r^2 / ( beta^2 + r^2 ),      r = t_s / t_n.
//
//  Cas limites : traction pure (t_s = 0) -> m = 0 -> G_C = G_Ic ;
//  cisaillement pur (t_n <= 0) -> m = 1 -> G_C = G_IIc.
// ---------------------------------------------------------------------
inline double bkMixFromTraction(double tn, double ts, double beta)
{
    double a = mac(tn);
    if (!(a > 0.0)) return 1.0;                 // cisaillement pur
    double r = std::fabs(ts) / a;
    double d = beta * beta + r * r;
    return (d > 0.0) ? (r * r) / d : 0.0;
}

// La forme litterale de l'eq. 18, sur les separations. Fournie pour les
// bancs de verification (elle doit redonner bkMixFromTraction le long
// d'un trajet proportionnel), jamais appelee dans la boucle de calcul.
inline double bkMixFromSeparation(double dn, double ds, double beta)
{
    double a = mac(dn), b = beta * ds;
    double d = a * a + b * b;
    return (d > 0.0) ? (b * b) / d : 0.0;
}

// ---------------------------------------------------------------------
//  §2.4 eq. 18 — ENERGIE DE RUPTURE MIXTE DE BENZEGGAGH-KENANE
//
//      G_C = G_Ic + (G_IIc - G_Ic) * m^eta,    eta dans [1,5 ; 2,5]
//
//  Evaluee au ratio de mode A L'INSERTION puis GELEE (« … puis gele »).
//  Le gel est essentiel : reevaluer G_C a chaque pas pendant
//  l'adoucissement fait varier delta_m^f sous les pieds de la variable
//  d'endommagement D = delta_m^max / delta_m^f, et l'integrale de la loi
//  cesse de valoir G_C.
// ---------------------------------------------------------------------
inline double bkEnergy(double GIc, double GIIc, double m, double eta)
{
    double mm = std::min(1.0, std::max(0.0, m));
    double e  = (eta > 0.0) ? eta : 1.0;
    return GIc + (GIIc - GIc) * std::pow(mm, e);
}

// ---------------------------------------------------------------------
//  ETAT FIGE A L'INSERTION
//
//  Tout ce que la loi initialement rigide doit retenir de l'instant
//  d'activation. Volontairement minimal : quatre doubles par joint.
// ---------------------------------------------------------------------
struct Stamp {
    double tmIns = 0.0;   // traction effective REELLEMENT transmise  [Pa]
    double beta  = 1.0;   // (t_s0 / t_n0) a l'insertion               [-]
    double Gc    = 0.0;   // G_C de Benzeggagh-Kenane, gele       [J/m^2]
    double dmF   = 0.0;   // longueur de rupture = 2 G_C / t_m^ins     [m]
    // ---- BRANCHE ASCENDANTE COURTE (§2.4, « Discontinuite temporelle des
    // lois extrinseques », Papoulia, Sam & Vavasis 2003) ---------------
    // dm0 = rise * dmF est un DECALAGE de separation effective, tamponne a
    // l'insertion : le joint nait AU SOMMET d'une branche elastique fictive
    // de longueur dm0, et non a l'origine. Trois consequences, toutes
    // necessaires — le run du 2026-09-11 (out_note2026, _dtsafe, _fricoff)
    // a prouve qu'a rise = 0 la loi est INCONDITIONNELLEMENT instable :
    //  (1) la raideur de charge/decharge est BORNEE a k0 = t_ins / dm0
    //      = t_ins^2 / (2 rise G_C). A rise = 0 elle est t_ins / dm_max,
    //      INFINIE pour un joint qui s'est ouvert d'un rien puis recharge :
    //      mesure a la trame 9 de out_note2026_dtsafe, p99,9 = 210 x pj et
    //      max = 549 x pj — hors de tout budget CFL, quel que soit dt (la
    //      division de dt par 3,49 n'a fait que RETARDER l'explosion) ;
    //  (2) aucun saut de traction a l'insertion : a separation geometrique
    //      nulle, la separation effective vaut dm0 dans la DIRECTION de la
    //      traction tamponnee (en, es), et l'eq. 19 rend exactement
    //      (t_n^ins, t_s^ins). A rise = 0, split() rendait ZERO a
    //      l'insertion — un Dirac de -t_ins sur chaque facette activee ;
    //  (3) la compliance ajoutee est negligeable (la note : « delta_m^0 =
    //      1e-3 delta_m^f »), et int t_m d(delta_m) sur la branche
    //      ADOUCISSANTE reste exactement G_C, car l'eq. 17 est ecrite dans
    //      la separation GEOMETRIQUE ; seule la decharge lit le decalage.
    double dm0   = 0.0;   // decalage = rise * dmF                       [m]
    double en    = 1.0;   // direction de la traction tamponnee dans
    double es    = 0.0;   // l'espace effectif (<t_n>, t_s/beta) / t_m  [-]
    bool   ok() const { return tmIns > 0.0 && dmF > 0.0; }
    // Raideur de charge de la branche ascendante, celle qui doit entrer au
    // budget CFL. 0 signifie « NON BORNEE » (rise = 0) : l'appelant doit
    // refuser ou avertir, jamais budgeter 0.
    double loadingStiffness() const {
        return (dm0 > 0.0 && tmIns > 0.0) ? tmIns / dm0 : 0.0;
    }
};

// Estimation A PRIORI de cette raideur pour une facette encore LIEE (elle
// peut s'inserer a tout pas, avec t_ins <= seuil x DIF) : c'est ce que le
// budget CFL doit porter pour TOUTES les facettes, liees comprises. Le run
// du 2026-09-11 a montre que les exclure du budget (schema du lot B) explose.
//   k0 = ft^2 / (2 rise Gf)
inline double loadingStiffnessEstimate(double ft, double Gf, double rise)
{
    return (ft > 0.0 && Gf > 0.0 && rise > 0.0) ? ft * ft / (2.0 * rise * Gf) : 0.0;
}

// ---------------------------------------------------------------------
//  §2.4 — TAMPON D'INSERTION
//
//  t_m^ins = t_m(sigma_F) : la traction REELLEMENT TRANSMISE par la
//  facette au pas d'activation, et NON le seuil nominal. C'est le point
//  (iii) de la section 2.4, celui qui supprime tout saut de contrainte,
//  « depassement explicite compris » : le joint nait exactement sur la
//  traction que la facette liee portait au pas precedent.
//
//      delta_m^f = 2 G_C / t_m^ins        (eq. 17)
//
//  Par construction l'integrale de t_m d(delta_m) de 0 a delta_m^f vaut
//  1/2 * t_m^ins * delta_m^f = G_C, INDEPENDAMMENT du maillage et du
//  depassement a l'insertion. C'est la propriete qui rend la loi
//  objective la ou une loi ancree sur le seuil nominal ne l'est pas.
//
//  tn, ts    : tractions transmises par la facette a l'activation [Pa]
//  tn0, ts0  : seuils dynamiques de la facette (DIF deja applique)  [Pa]
//  GIc, GIIc : energies de rupture de la facette                [J/m^2]
//  eta       : exposant de Benzeggagh-Kenane ; eta <= 0 => pas de
//              melange, G_C = G_Ic (jointMixLaw = none)
// ---------------------------------------------------------------------
//  rise      : longueur de la branche ascendante en FRACTION de delta_m^f
//              (§2.4 : 1e-3). 0 = pas de branche = loi litteralement rigide,
//              INSTABLE en explicite (voir Stamp) ; a n'utiliser que sur un
//              banc qui veut justement le montrer.
// ---------------------------------------------------------------------------
inline Stamp stampInsertion(double tn, double ts, double tn0, double ts0,
                            double GIc, double GIIc, double eta,
                            double rise = 0.0)
{
    Stamp S;
    S.beta = (tn0 > 0.0 && ts0 > 0.0) ? ts0 / tn0 : 1.0;
    if (!(S.beta > 0.0)) S.beta = 1.0;
    S.tmIns = effTrac(tn, ts, S.beta);
    if (eta > 0.0) {
        double m = bkMixFromTraction(tn, ts, S.beta);
        S.Gc = bkEnergy(GIc, GIIc, m, eta);
    } else {
        S.Gc = GIc;
    }
    S.dmF = (S.tmIns > 0.0) ? 2.0 * S.Gc / S.tmIns : 0.0;
    // direction de la traction tamponnee dans l'espace effectif
    // (<t_n>, t_s/beta) / t_m : c'est la direction dans laquelle le
    // decalage dm0 est pose, pour que l'eq. 19 rende (t_n^ins, t_s^ins)
    // exactement a separation geometrique nulle
    if (S.tmIns > 0.0) {
        S.en = mac(tn) / S.tmIns;
        S.es = (std::fabs(ts) / S.beta) / S.tmIns;
    }
    S.dm0 = (rise > 0.0) ? rise * S.dmF : 0.0;
    return S;
}

// ---------------------------------------------------------------------------
//  SEPARATION EFFECTIVE AVEC DECALAGE
//
//  Entree : la separation GEOMETRIQUE (dn, ds) — ds est la NORME du
//  glissement, la direction tangentielle est geree par l'appelant.
//  Sortie : la separation effective decalee, dans l'espace (dn, beta ds) :
//      dnEff        = <dn> + dm0 en
//      (beta ds)Eff = beta ds + dm0 es          ->  dsEff = ds + dm0 es / beta
//  et sa norme dmEff = sqrt(dnEff^2 + beta^2 dsEff^2), qui vaut EXACTEMENT
//  dm0 a separation geometrique nulle (en^2 + es^2 = 1).
//  A rise = 0 : dnEff = <dn>, dsEff = ds, dmEff = effSep(dn, ds, beta).
//
//  ATTENTION a la compression : <dn> ecrete la partie normale. La
//  compression est reprise par la penalite de contact du solveur (§2.4),
//  jamais par la loi cohesive ; le decalage ne s'y applique pas.
// ---------------------------------------------------------------------------
inline void effOffset(const Stamp& S, double dn, double ds,
                      double& dnEff, double& dsEff, double& dmEff)
{
    dnEff = mac(dn) + S.dm0 * S.en;
    dsEff = ds + (S.beta > 0.0 ? S.dm0 * S.es / S.beta : 0.0);
    double b = S.beta * dsEff;
    dmEff = std::sqrt(dnEff * dnEff + b * b);
}

// ---------------------------------------------------------------------
//  §2.4 eq. 17 — LOI INITIALEMENT RIGIDE, DECHARGE SECANTE
//
//      charge    : t_m = t_m^ins ( 1 - delta_m / delta_m^f )
//      decharge  : t_m = t_m^ins ( 1 - D ) delta_m / delta_m^max
//      D         = delta_m^max / delta_m^f   dans [0, 1]
//
//  AUCUNE raideur initiale K_n : la branche elastique pre-rupture est
//  portee par la MATRICE seule (point (i) de la section 2.4). C'est de
//  la que vient le gain sur le pas de temps : la penalite des facettes
//  encore liees ne charge plus le budget CFL.
//
//  dm    : separation effective courante             [m]
//  dmMax : plus grande separation effective atteinte [m]
//  Retourne la traction effective t_m >= 0 [Pa].
// ---------------------------------------------------------------------
//  Les deux arguments sont des separations EFFECTIVES (sorties de effOffset,
//  ou la separation geometrique a rise = 0). L'enveloppe de l'eq. 17 est
//  ecrite dans la separation GEOMETRIQUE dm - dm0, de sorte que l'integrale
//  sur la branche adoucissante vaut exactement G_C quel que soit rise ; la
//  decharge, elle, est la secante vers l'ORIGINE EFFECTIVE, de pente bornee
//  par t_ins / dm0. A rise = 0 (dm0 = 0) tout se reduit a la forme initiale.
// ---------------------------------------------------------------------------
inline double traction(const Stamp& S, double dmEff, double dmMaxEff)
{
    if (!S.ok()) return 0.0;
    double dmx = std::max(dmMaxEff, dmEff);
    double g   = dmx - S.dm0;                            // ouverture geometrique max
    if (g >= S.dmF) return 0.0;                          // joint rompu
    double env = S.tmIns * (1.0 - mac(g) / S.dmF);       // enveloppe au pic
    if (dmEff >= dmx) return env;                        // en charge
    return (dmx > 0.0) ? env * (dmEff / dmx) : 0.0;      // secante a l'origine
}

// Variable d'endommagement de l'eq. 17, D = delta_m^max / delta_m^f, sur
// l'ouverture GEOMETRIQUE (le decalage n'endommage pas).
inline double damage(const Stamp& S, double dmMaxEff)
{
    if (!S.ok()) return 0.0;
    return std::min(1.0, std::max(0.0, (dmMaxEff - S.dm0) / S.dmF));
}

// ---------------------------------------------------------------------
//  §2.4 eq. 19 — PARTITION DE LA TRACTION EFFECTIVE
//
//      t_n       = (t_m / delta_m) <delta_n>
//      t_s,alpha = (t_m / delta_m) beta^2 delta_s,alpha
//
//  `tn` recoit la part normale (0 en compression : la compression est
//  reprise par le contact penalise, cf. §2.4 « compression et
//  frottement sur un joint insere »), `tsScale` le FACTEUR a appliquer
//  au vecteur de glissement pour obtenir la traction tangentielle.
// ---------------------------------------------------------------------
inline void split(double tm, double dm, double dn, double beta,
                  double& tn, double& tsScale)
{
    if (!(dm > 0.0) || !(tm > 0.0)) { tn = 0.0; tsScale = 0.0; return; }
    double k = tm / dm;
    tn      = k * mac(dn);
    tsScale = k * beta * beta;
}

// ---------------------------------------------------------------------
//  §2.4 — MELANGE COHESIF / FROTTANT SUR UN JOINT EN COMPRESSION
//
//      tau_lim = (1 - D) * c_coh  +  D * mu * <-t_n>
//
//  Le facteur D sur la part FROTTANTE est le point que rockim_g0 n'avait
//  pas : sans lui le frottement vaut mu<-t_n> a pleine valeur des D = 0,
//  et le cisaillement est compte deux fois — une fois par la
//  viscoplasticite de la matrice, une fois par le joint. C'est
//  exactement ce que la note interdit (« un mecanisme par physique »).
//
//  `mobilised = false` restitue le comportement historique (frottement
//  plein des D = 0), pour que le defaut reste bit-identique.
// ---------------------------------------------------------------------
//  Deux ENVELOPPES pour le terme frottant, selon la cle jointEnvelope :
//    yan  (defaut) : mu <-t_n>              — Yan et al. 2023
//    yang          : mu (-min(t_n, ft))     — Yang et al. 2025, leur eq. 1,
//                    qui garde une branche en traction jusqu'a ft.
//  Le facteur D est le meme dans les deux cas. Les deux solveurs
//  appelaient mcFrictionTerm() a la main pour le cas yang : c'est la meme
//  expression, ramenee ici pour qu'il n'en existe qu'une.
// ---------------------------------------------------------------------------
inline double shearCap(double cohEff, double mu, double tnNeg, double D,
                       bool mobilised, bool yangEnv = false, double ft = 0.0)
{
    double term = yangEnv ? -std::min(tnNeg, ft) : mac(-tnNeg);
    double fric = mu * term;
    double lim  = cohEff + (mobilised ? D * fric : fric);
    return lim > 0.0 ? lim : 0.0;
}

// ---------------------------------------------------------------------
//  §2.5 — OPTION B : TERME VISQUEUX APRES INSERTION
//
//      t = t_coh(delta_m) + eta * d(delta)/dt
//
//  applicable SEULEMENT apres insertion. La note insiste : UNE SEULE des
//  deux options (DIF a l'insertion OU viscosite apres insertion), jamais
//  les deux — l'exclusion est verifiee a la lecture des cles par le
//  solveur, pas ici.
// ---------------------------------------------------------------------
inline double viscous(double eta, double rate) { return eta * rate; }

// ---------------------------------------------------------------------
//  §1.3 — LONGUEUR DE BANDE DE FISSURATION D'UN TETRAEDRE
//
//      h_e = ( 12 V_e / sqrt(2) )^(1/3)
//
//  C'est l'arete d'un tetraedre REGULIER de volume V_e. rockim_g0
//  passait V^(1/3) en fem3d (facteur 0,49 vs la note) et le diametre
//  inscrit 6V/A en fdem3d (facteur 0,20) : un G_c^bulk identifie avec la
//  convention de la note y produisait une branche adoucissante 2 a 5
//  fois trop ductile. Pose ici pour que les trois solveurs partagent la
//  meme definition.
// ---------------------------------------------------------------------
inline double tetEdgeLength(double V)
{
    return (V > 0.0) ? std::cbrt(12.0 * V / std::sqrt(2.0)) : 0.0;
}

// ---------------------------------------------------------------------
//  §1.3 eq. 7 — CONDITION DE VALIDITE DE LA BANDE DE FISSURATION
//
//      h_e  <  2 E G_c^bulk / sigma_c0^2
//
//  (pas de snap-back elementaire ; ~22 mm pour Kuru avec
//  G_c^bulk = 10 kJ/m^2). Retourne la borne ; l'appelant decide s'il
//  leve une exception ou un avertissement.
// ---------------------------------------------------------------------
inline double compBandLimit(double E, double GcBulk, double sigC0)
{
    return (sigC0 > 0.0) ? 2.0 * E * GcBulk / (sigC0 * sigC0) : 0.0;
}

}  // namespace jtsl
}  // namespace rockim
