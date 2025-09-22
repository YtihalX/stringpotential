#pion decay constant
const f_pi = 0.092 #GeV
#and coupling
const g_pi = 0.5704 #*2.0 #GeV

#Goldstone boson masses
const m_pi = 0.138039407 #GeV
const m_K = 0.498 #GeV
const m_eta = 0.548 #GeV

#heavy meson masses
const m_B = 5.27934 #GeV
const m_B_star = 5.32471 #GeV
const m_B_s = 5.36692 #GeV
const m_B_star_s = 5.4154 #GeV

#asign masses to channels
const m11 = m_B
const m12 = m_B_star

const m21 = m_B_s
const m22 = m_B_star_s
const m = [m11 m12; m21 m22]

#masses of bound states in 1++ channel
const m_Xb14P = 10.825 - (m11 + m12)
const m_Xb13P = 10.529 - (m11 + m12)
const m_Xb12P = 10.224 - (m11 + m12)
const m_Xb11P = 9.909 - (m11 + m12)

#use couplings from D200 fit
#Lattice spacing
const a_l = 0.06426 * 5.068 #GeV^(-1)
#pion mass dependet scale factor
const g = 0.898794378386677 / a_l#*1e1 #increase coupling just for illustration
#channel dependent couplings
const g1 = g * 0.014926616931653945
const g2 = g * 0.006467550544943349
const delta = [0, 0, m21 + m22 - m11 - m12]
const mu = [m11 * m12 / (m11 + m12), m11 * m12 / (m11 + m12), m21 * m22 / (m21 + m22)]
const ϵ = 1e-9
#constants to express Schroedinger equation in dimensionless units (according to Cornell paper)
m_c = 1.85 #GeV
a_Cornell = 1.95 #GeV^-1
