use std::{collections::HashMap};

// Hashmap of all element and their atomic number
pub fn elem_to_idx(elem: &str) -> usize {
    let elem_list:HashMap<String, usize> = HashMap::from([
        (String::from("H"),1),
        (String::from("He"),2),
        (String::from("Li"),3),
        (String::from("Be"),4),
        (String::from("B"),5),
        (String::from("C"),6),
        (String::from("N"),7),
        (String::from("O"),8),
        (String::from("F"),9),
        (String::from("Ne"),10),
        (String::from("Na"),11),
        (String::from("Mg"),12),
        (String::from("Al"),13),
        (String::from("Si"),14),
        (String::from("P"),15),
        (String::from("S"),16),
        (String::from("Cl"),17),
        (String::from("Ar"),18),
        (String::from("K"),19),
        (String::from("Ca"),20),
        (String::from("Sc"),21),
        (String::from("Ti"),22),
        (String::from("V"),23),
        (String::from("Cr"),24),
        (String::from("Mn"),25),
        (String::from("Fe"),26),
        (String::from("Co"),27),
        (String::from("Ni"),28),
        (String::from("Cu"),29),
        (String::from("Zn"),30),
        (String::from("Ga"),31),
        (String::from("Ge"),32),
        (String::from("As"),33),
        (String::from("Se"),34),
        (String::from("Br"),35),
        (String::from("Kr"),36),
        (String::from("Rb"),37),
        (String::from("Sr"),38),
        (String::from("Y"),39),
        (String::from("Zr"),40),
        (String::from("Nb"),41),
        (String::from("Mo"),42),
        (String::from("Tc"),43),
        (String::from("Ru"),44),
        (String::from("Rh"),45),
        (String::from("Pd"),46),
        (String::from("Ag"),47),
        (String::from("Cd"),48),
        (String::from("In"),49),
        (String::from("Sn"),50),
        (String::from("Sb"),51),
        (String::from("Te"),52),
        (String::from("I"),53),
        (String::from("Xe"),54),
        (String::from("Cs"),55),
        (String::from("Ba"),56),
        (String::from("La"),57),
        (String::from("Ce"),58),
        (String::from("Pr"),59),
        (String::from("Nd"),60),
        (String::from("Pm"),61),
        (String::from("Sm"),62),
        (String::from("Eu"),63),
        (String::from("Gd"),64),
        (String::from("Tb"),65),
        (String::from("Dy"),66),
        (String::from("Ho"),67),
        (String::from("Er"),68),
        (String::from("Tm"),69),
        (String::from("Yb"),70),
        (String::from("Lu"),71),
        (String::from("Hf"),72),
        (String::from("Ta"),73),
        (String::from("W"),74),
        (String::from("Re"),75),
        (String::from("Os"),76),
        (String::from("Ir"),77),
        (String::from("Pt"),78),
        (String::from("Au"),79),
        (String::from("Hg"),80),
        (String::from("Tl"),81),
        (String::from("Pb"),82),
        (String::from("Bi"),83),
        (String::from("Po"),84),
        (String::from("At"),85),
        (String::from("Rn"),86),
        (String::from("Fr"),87),
        (String::from("Ra"),88),
        (String::from("Ac"),89),
        (String::from("Th"),90),
        (String::from("Pa"),91),
        (String::from("U"),92),
        (String::from("Np"),93),
        (String::from("Pu"),94),
        (String::from("Am"),95),
        (String::from("Cm"),96),
        (String::from("Bk"),97),
        (String::from("Cf"),98),
        (String::from("Es"),99),
        (String::from("Fm"),100),
        (String::from("Md"),101),
        (String::from("No"),102),
        (String::from("Lr"),103),
        (String::from("Rf"),104),
        (String::from("Db"),105),
        (String::from("Sg"),106),
        (String::from("Bh"),107),
        (String::from("Hs"),108),
        (String::from("Mt"),109),
        (String::from("Ds"),110),
        (String::from("Rg"),111),
        (String::from("Cn"),112),
        (String::from("Nh"),113),
        (String::from("Fl"),114),
        (String::from("Mc"),115),
        (String::from("Lv"),116),
        (String::from("Ts"),117),
        (String::from("Og"),118)]);

        elem_list.get(elem).unwrap().clone()
    
}

pub fn elems_to_idx(elems: &Vec<String>) -> Vec<usize> {
    let mut idxs = Vec::new();
    for elem in elems {
        idxs.push(elem_to_idx(elem));
    }
    idxs
}

// [charge, mass, name]
const ELEM_INFO: [[f64; 2]; 119] = [
    [0.0, 0.0],         // 0
    [1.0, 1.00784],     // 1
    [2.0, 4.002602],   // 2
    [3.0, 6.938],      // 3
    [4.0, 9.0121831],  // 4
    [5.0, 10.81],       // 5
    [6.0, 12.011],      // 6
    [7.0, 14.007],      // 7
    [8.0, 15.999],      // 8
    [9.0, 18.998],      // 9
    [10.0, 20.180],    // 10
    [11.0, 22.990],    // 11
    [12.0, 24.305],    // 12
    [13.0, 26.982],    // 13
    [14.0, 28.085],    // 14
    [15.0, 30.974],     // 15
    [16.0, 32.06],      // 16
    [17.0, 35.45],     // 17
    [18.0, 39.948],    // 18
    [19.0, 39.098],     // 19
    [20.0, 40.078],    // 20
    [21.0, 44.956],    // 21
    [22.0, 47.867],    // 22
    [23.0, 50.942],     // 23
    [24.0, 51.996],    // 24
    [25.0, 54.938],    // 25
    [26.0, 55.845],    // 26
    [27.0, 58.933],    // 27
    [28.0, 58.693],    // 28
    [29.0, 63.546],    // 29
    [30.0, 65.38],     // 30
    [31.0, 69.723],    // 31
    [32.0, 72.63],     // 32
    [33.0, 74.922],    // 33
    [34.0, 78.971],    // 34
    [35.0, 79.904],    // 35
    [36.0, 83.798],    // 36
    [37.0, 85.468],    // 37
    [38.0, 87.62],     // 38
    [39.0, 88.906],     // 39
    [40.0, 91.224],    // 40
    [41.0, 92.906],    // 41
    [42.0, 95.95],     // 42
    [43.0, 98.907],    // 43
    [44.0, 101.07],    // 44
    [45.0, 102.91],    // 45
    [46.0, 106.42],    // 46
    [47.0, 107.87],    // 47
    [48.0, 112.41],    // 48
    [49.0, 114.82],    // 49
    [50.0, 118.71],    // 50
    [51.0, 121.76],    // 51
    [52.0, 127.6],     // 52
    [53.0, 126.9],      // 53
    [54.0, 131.29],    // 54
    [55.0, 132.91],    // 55
    [56.0, 137.33],    // 56
    [57.0, 138.91],    // 57
    [58.0, 140.12],    // 58
    [59.0, 140.91],    // 59
    [60.0, 144.24],    // 60
    [61.0, 145.0],     // 61
    [62.0, 150.36],    // 62
    [63.0, 151.96],    // 63
    [64.0, 157.25],    // 64
    [65.0, 158.93],    // 65
    [66.0, 162.5],     // 66
    [67.0, 164.93],    // 67
    [68.0, 167.26],    // 68
    [69.0, 168.93],    // 69
    [70.0, 173.04],    // 70
    [71.0, 174.97],    // 71
    [72.0, 178.49],    // 72
    [73.0, 180.95],    // 73
    [74.0, 183.84],     // 74
    [75.0, 186.21],    // 75
    [76.0, 190.23],    // 76
    [77.0, 192.22],    // 77
    [78.0, 195.08],    // 78
    [79.0, 196.97],    // 79
    [80.0, 200.59],    // 80
    [81.0, 204.38],    // 81
    [82.0, 207.2],     // 82
    [83.0, 208.98],    // 83
    [84.0, 208.98],    // 84
    [85.0, 209.99],    // 85
    [86.0, 222.02],    // 86
    [87.0, 223.02],    // 87
    [88.0, 226.03],    // 88
    [89.0, 227.03],    // 89
    [90.0, 232.04],    // 90
    [91.0, 231.04],    // 91
    [92.0, 238.03],     // 92
    [93.0, 237.05],    // 93
    [94.0, 244.06],    // 94
    [95.0, 243.06],    // 95
    [96.0, 247.07],    // 96
    [97.0, 247.07],    // 97
    [98.0, 251.08],    // 98
    [99.0, 252.08],    // 99
    [100.0, 257.1],    // 100
    [101.0, 258.1],    // 101
    [102.0, 259.1],    // 102
    [103.0, 262.11],   // 103
    [104.0, 267.12],   // 104
    [105.0, 270.13],   // 105
    [106.0, 271.13],   // 106
    [107.0, 270.13],   // 107
    [108.0, 277.15],   // 108
    [109.0, 278.15],   // 109
    [110.0, 281.16],   // 110
    [111.0, 282.16],   // 111
    [112.0, 285.17],   // 112
    [113.0, 286.17],   // 113
    [114.0, 289.18],   // 114
    [115.0, 290.18],   // 115
    [116.0, 293.19],   // 116
    [117.0, 294.19],   // 117
    [118.0, 294.19]];   // 118
    
    
pub const R0_CC: f64 = 1.5300;
pub const R0_CH: f64 = 1.1100;
pub const Kb_CC: f64 = 300.0;
pub const Kb_CH: f64 = 350.0;
pub const A0_CCC: f64 = 109.50;
pub const A0_XCX: f64 = 109.50;
pub const Ka_CCC: f64 = 60.00;
pub const Ka_XCX: f64 = 35.00;
pub const N_XCCX: f64 = 3.0;
pub const A_XCCX: f64 = 0.300;
pub const EPSILON_H: f64 = 0.0300;
pub const EPSILON_C: f64 = 0.0700;
pub const SIGMA_H: f64 = 1.20000;
pub const SIGMA_C: f64 = 1.75000;



