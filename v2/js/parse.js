const window_size = 25;
const amino_acids = {
    A: { name: 'ala', ppii_h: 0.37, helix: 1.42, hydr: 0.0728, lambda_1: 0.0, lambda_2: 0.0 },
    C: { name: 'cys', ppii_h: 0.25, helix: 0.73, hydr: 0.3557, lambda_1: 0.0, lambda_2: 0.0 },
    D: { name: 'asp', ppii_h: 0.3, helix: 1.01, hydr: -0.0552, lambda_1: -1.0, lambda_2: -1.0 },
    E: { name: 'glu', ppii_h: 0.42, helix: 1.63, hydr: -0.0295, lambda_1: -1.0, lambda_2: -1.0 },
    F: { name: 'phe', ppii_h: 0.17, helix: 1.16, hydr: 0.4201, lambda_1: 0.0, lambda_2: 0.0 },
    G: { name: 'gly', ppii_h: 0.13, helix: 0.50, hydr: -0.0589, lambda_1: 0.0, lambda_2: 0.0 },
    H: { name: 'his', ppii_h: 0.2, helix: 1.20, hydr: 0.0874, lambda_1: 0.0, lambda_2: 0.0 },
    I: { name: 'ile', ppii_h: 0.39, helix: 1.12, hydr: 0.3805, lambda_1: 0.0, lambda_2: 0.0 },
    K: { name: 'lys', ppii_h: 0.56, helix: 1.24, hydr: -0.0053, lambda_1: 1.0, lambda_2: 1.0 },
    L: { name: 'leu', ppii_h: 0.24, helix: 1.29, hydr: 0.3819, lambda_1: 0.0, lambda_2: 0.0 },
    M: { name: 'met', ppii_h: 0.36, helix: 1.21, hydr: 0.1613, lambda_1: 0.0, lambda_2: 0.0 },
    N: { name: 'asn', ppii_h: 0.27, helix: 0.71, hydr: -0.0390, lambda_1: 0.0, lambda_2: 0.0 },
    P: { name: 'pro', ppii_h: 1.0, helix: 0.65, hydr: -0.0492, lambda_1: 0.0, lambda_2: 0.0 },
    Q: { name: 'gln', ppii_h: 0.53, helix: 1.02, hydr: 0.0126, lambda_1: 0.0, lambda_2: 0.0 },
    R: { name: 'arg', ppii_h: 0.38, helix: 1.06, hydr: 0.0394, lambda_1: 1.0, lambda_2: 1.0 },
    S: { name: 'ser', ppii_h: 0.24, helix: 0.71, hydr: -0.0282, lambda_1: 0.0, lambda_2: 0.0 },
    T: { name: 'thr', ppii_h: 0.32, helix: 0.78, hydr: 0.0239, lambda_1: 0.0, lambda_2: 0.0 },
    V: { name: 'val', ppii_h: 0.39, helix: 0.99, hydr: 0.2947, lambda_1: 0.0, lambda_2: 0.0 },
    W: { name: 'trp', ppii_h: 0.25, helix: 1.05, hydr: 0.4114, lambda_1: 0.0, lambda_2: 0.0 },
    Y: { name: 'tyr', ppii_h: 0.25, helix: 0.67, hydr: 0.3113, lambda_1: 0.0, lambda_2: 0.0 },
};
const amino_acid_keys = Object.keys(amino_acids);

function getResidues(input) {
    var residues = [];
    for (let i = 0; i < input.length; i++) {
        residues.push({
            amino_acid: input[i],
            region: null,
            region_pi_q: null,
            dist_norm: null,
            dist_norm_pi_q: null,
        });
    }

    for (let i = 0; i < input.length; i++) {
        if (i < input.length - window_size + 1) {
            var middle_position = i + Math.trunc(window_size / 2);
            var count = {
                A: 0,
                C: 0,
                D: 0,
                E: 0,
                F: 0,
                G: 0,
                H: 0,
                I: 0,
                K: 0,
                L: 0,
                M: 0,
                N: 0,
                P: 0,
                Q: 0,
                R: 0,
                S: 0,
                T: 0,
                V: 0,
                W: 0,
                Y: 0,
            };

            for (let j = i; j < i + window_size; j++) {
                count[input[j]]++;
            }

            var ppii = 0.0;
            var helix = 0.0;
            var hydr = 0.0;

            amino_acid_keys.forEach((x) => {
                ppii += count[x] * amino_acids[x].ppii_h;
                helix += count[x] * amino_acids[x].helix;
                hydr += count[x] * amino_acids[x].hydr;
            });
            ppii = ppii / window_size;
            helix = helix / window_size;
            hydr = hydr / window_size;

            if (hydr >= 0.08280152) {
                residues[middle_position].region = 'F';
                residues[middle_position].region_pi_q = residues[middle_position].region;
                residues[middle_position].dist_norm = (hydr - 0.08280152) / (0.01679414 * 2.0);
                residues[middle_position].dist_norm_pi_q = residues[middle_position].dist_norm;
            } else {
                var RY = count.R * count.Y / (count.R == count.Y ? 1 : Math.abs(count.R - count.Y));
                var RF = count.R * count.F / (count.R == count.F ? 1 : Math.abs(count.R - count.F));
                var KY = count.K * count.Y / (count.K == count.Y ? 1 : Math.abs(count.K - count.Y));
                var KF = count.K * count.F / (count.K == count.F ? 1 : Math.abs(count.K - count.F));
                var FY = count.F * count.Y / (count.F == count.Y ? 1 : Math.abs(count.F - count.Y));

                var U_pi = 0.137 * (3.0 * RY + 2.0 * KY + 2.0 * RF + 1.0 * KF + 1.0 * FY);

                var scd = 0;

                for (let j = i; j < i + window_size; j++) {
                    for (let k = i; k < i + window_size; k++) {
                        if (k > j) {

                            var lambda_1 = amino_acids[input[j]].lambda_1;
                            var lambda_2 = amino_acids[input[k]].lambda_2;

                            scd += (lambda_1 * lambda_2) * Math.sqrt(Math.abs(k - j));
                        }
                    }
                }

                scd = scd / window_size

                var net_charge = Math.abs(count.D + count.E - count.K - count.R);
                var ncpr = net_charge / window_size

                var U_q = 8.4 * scd + 5.6 * ncpr;

                if (ppii == 1.0) {
                    ppii = 0.98;
                }

                var v_exponent = 0.503 - 0.11 * Math.log(1.0 - ppii);
                var rh = 2.16 * (4 * window_size) ** v_exponent + 0.26 * 4 * net_charge - 0.29 * (4 * window_size) ** 0.5;
                var nu_model = Math.log(rh / 2.16) / Math.log(4 * window_size);

                var c1 = -0.244078945;
                var c2 = 0.7885823;
                var c3 = 0.5582901;
                var c4 = 1.022552;

                var m = -1.0 / c1
                var b = c3 - m * c4
                var x = (b - c2) / (c1 - m)
                var y = m * x + b
                var id_dist = Math.sqrt((c4 - x) * (c4 - x) + (c3 - y) * (c3 - y))

                c3 = 0.5416;
                c4 = 0.9327272;

                b = c3 - m * c4
                x = (b - c2) / (c1 - m)
                y = m * x + b
                var ps_dist = Math.sqrt((c4 - x) * (c4 - x) + (c3 - y) * (c3 - y))

                b = nu_model - m * helix
                x = (b - c2) / (c1 - m)
                y = m * x + b

                if (((nu_model - c2) / c1) <= helix) {
                    residues[middle_position].region = 'D';
                    residues[middle_position].region_pi_q = residues[middle_position].region;
                    residues[middle_position].dist_norm = Math.sqrt((helix - x) * (helix - x) + (nu_model - y) * (nu_model - y)) / id_dist;
                    residues[middle_position].dist_norm_pi_q = residues[middle_position].dist_norm;

                    if (residues[middle_position].dist_norm < U_pi + U_q) {
                        residues[middle_position].region_pi_q = 'P';
                        residues[middle_position].dist_norm_pi_q = U_pi + U_q - residues[middle_position].dist_norm;
                    }
                } else {
                    residues[middle_position].region = 'P';
                    residues[middle_position].region_pi_q = residues[middle_position].region;
                    residues[middle_position].dist_norm = Math.sqrt((helix - x) * (helix - x) + (nu_model - y) * (nu_model - y)) / ps_dist;
                    residues[middle_position].dist_norm_pi_q = U_pi + U_q + residues[middle_position].dist_norm;
                }
            }
        }
    }

    for (let i = 0; i < Math.trunc(window_size / 2); i++) {
        residues[i].region = residues[Math.trunc(window_size / 2)].region;
        residues[i].region_pi_q = residues[Math.trunc(window_size / 2)].region_pi_q;
        residues[i].dist_norm = residues[Math.trunc(window_size / 2)].dist_norm;
        residues[i].dist_norm_pi_q = residues[Math.trunc(window_size / 2)].dist_norm_pi_q;

        residues[input.length - i - 1].region = residues[input.length - Math.trunc(window_size / 2) - 1].region;
        residues[input.length - i - 1].region_pi_q = residues[input.length - Math.trunc(window_size / 2) - 1].region_pi_q;
        residues[input.length - i - 1].dist_norm = residues[input.length - Math.trunc(window_size / 2) - 1].dist_norm;
        residues[input.length - i - 1].dist_norm_pi_q = residues[input.length - Math.trunc(window_size / 2) - 1].dist_norm_pi_q;
    }

    return residues;
}

function getDomainRegion(residues, regionLength, cutoff, region) {
    var i = 0;
    var regions = [];

    while (i + regionLength <= residues.length) {
        var j = i;
        var found_region = false;
        var count_p = 0;
        var count_w = 0;

        while (j < residues.length) {
            count_w++;
            if (residues[j].region == region) {
                count_p++;
            }

            if (count_w >= regionLength) {
                var percent_p = count_p / count_w;
                if (percent_p >= cutoff) {
                    found_region = true;
                } else {
                    break;
                }
            }

            j++;
        }

        if (found_region) {
            regions.push({
                region: region,
                start: i,
                end: j - 1,
            });
            i = j;
        } else {
            i++;
        }
    }

    return regions;
}

function getRegions(residues, regionLength, cutoff) {
    var regions = [];
    var regionList;

    regionList = getDomainRegion(residues, regionLength, cutoff, 'P');
    for (const x of regionList) {
        regions.push(x);
    }
    regionList = getDomainRegion(residues, regionLength, cutoff, 'D');
    for (const x of regionList) {
        regions.push(x);
    }
    regionList = getDomainRegion(residues, regionLength, cutoff, 'F');
    for (const x of regionList) {
        regions.push(x);
    }

    regions.sort((a, b) => {
        return a.start - b.start;
    });

    if (regions.length > 1) {
        for (let i = 0; i < regions.length - 1; i++) {
            if (regions[i].end >= regions[i + 1].start) {
                regions[i].end = regions[i + 1].start + Math.round(((1.0 - cutoff) * regionLength) / 2.0);
                regions[i + 1].start = regions[i].end + 1;
            }
        }
    }

    return regions;
}

function parse(sequence, regionLength, cutoff) {
    var result = {
        error: false,
        message: null,
        data: null,
    };

    if (!sequence) {
        result.error = true;
        result.message = `No sequence specified`;
        return result;
    }

    sequence = sequence.replace(/\s+/g, '');
    sequence = sequence.toUpperCase();

    sequence = Array.from(sequence);

    if (sequence.length < 25) {
        result.error = true;
        result.message = `Sequence length should be at least 25 characters long`;
        return result;
    } else if (sequence.length > 10000) {
        result.error = true;
        result.message = `Sequence length should be at most 10000 characters long`;
        return result;
    }

    var residues = getResidues(sequence);
    var regions = getRegions(residues, regionLength, cutoff);

    var data = {
        sequence: sequence,
        residues: residues,
        regions: regions,
    };

    result.data = data;

    return result;
}
