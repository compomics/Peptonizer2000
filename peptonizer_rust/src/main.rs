use wasm_bindgen::prelude::*;
use serde::Deserialize;
use serde_json::{Value, Result};
use std::collections::HashMap;

extern crate wasm_bindgen;
extern crate serde_json;
extern crate serde;

pub fn perform_taxa_weighing(
    unipept_responses: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: i32,
    taxa_rank: String
) {

    let mut unipept_responses: Vec<UnipeptJson> = serde_json::from_str(&unipept_responses).unwrap();

    let pep_scores_map: HashMap<String, f32> = serde_json::from_str(&pep_scores).unwrap();
    let mut pep_scores: Vec<f32> = vec![0.0; pep_scores_map.len()];
    for i in 0..pep_scores_map.len() {
        pep_scores[i] = pep_scores_map[&unipept_responses[i].sequence];
    }

    let pep_psm_counts_map: HashMap<String, i32> = serde_json::from_str(&pep_psm_counts).unwrap();
    let mut pep_psm_counts: Vec<i32> = vec![0; pep_psm_counts_map.len()];
    for i in 0..pep_psm_counts_map.len() {
        pep_psm_counts[i] = pep_psm_counts_map[&unipept_responses[i].sequence];
    }

    normalize_unipept_responses(&mut unipept_responses, &taxa_rank);
}

#[derive(Debug, Deserialize)]
struct UnipeptJson {
    sequence: String,
    taxa: Vec<i32>,
}

fn normalize_unipept_responses(unipept_responses: &mut Vec<UnipeptJson>, taxa_rank: &str) {
    // Map all taxa onto the rank specified by the user
    for i in 0..unipept_responses.len() {
        unipept_responses[i].taxa = get_unique_lineage_at_specified_rank(&unipept_responses[i].taxa, taxa_rank);
    }
}

fn get_unique_lineage_at_specified_rank(tax_ids: &Vec<i32>, taxa_rank: &str) -> Vec<i32> {
    return vec!(1);
}



pub fn main() {

    let pep_scores: String = "{\"GIVKDTTGEPVIGANVVVK\":0.999,\"DTTGEPVIGANVVVK\":0.999,\"GTTTGTITDFDGNFQLSAK\":0.999,\"KNDMTGSVMAIKPDELSK\":0.999,\"NDMTGSVMAIKPDELSK\":0.999,\"NDMTGSVMAIKPDELSKGITTNAQDMLSGK\":0.999,\"GITTNAQDMLSGK\":0.999,\"TTNAQDMLSGK\":0.999,\"IAGVSVISNDGTPGGGAQIR\":0.999,\"SVISNDGTPGGGAQIR\":0.999,\"ISNDGTPGGGAQIR\":0.999,\"SNDGTPGGGAQIR\":0.999,\"IRGGSSLNASNDPLIVIDGLAIDNEGIK\":0.999,\"GGSSLNASNDPLIVIDGLAIDNEGIK\":0.999,\"GMANGLSMVNPADIETLTVLK\":0.999,\"GMANGLSMVNPADIETLTVLKDASATAIYGSR\":0.999,\"DASATAIYGSR\":0.999,\"ASNGVIIITTK\":0.999,\"ASNGVIIITTKK\":0.999,\"SNGVIIITTKK\":0.999,\"NGQAPSVSYNGSVSFSK\":0.999,\"RYDVLSGDEYR\":0.999,\"YDVLSGDEYR\":0.999,\"AYANQLWGDK\":0.999,\"AYANQLWGDKLPADLGTANTDWQDQIFR\":0.999,\"LPADLGTANTDWQDQIFR\":0.999,\"TAVSTDHHVSINGGFK\":0.999,\"AVSTDHHVSINGGFK\":0.999,\"VSTDHHVSINGGFK\":0.999,\"STDHHVSINGGFK\":0.999,\"VSLGYTDDNGIVK\":0.999,\"FTASVNLAPSFFEDHLK\":0.999,\"ASVNLAPSFFEDHLK\":0.999,\"YADTGAAIGGALAIDPTRPVYSNEDPYQFTGGYWQNINSTTGFSNPDWK\":0.999,\"YTSNPNSPQNPLAALELKNDKANSNDFVGNVDVDYK\":0.999,\"YTSNPNSPQNPLAALELK\":0.999,\"YTSNPNSPQNPLAALELKNDK\":0.999,\"PNSPQNPLAALELK\":0.999,\"SPQNPLAALELK\":0.999,\"NDKANSNDFVGNVDVDYK\":0.999,\"NDKANSNDFVGNVDVDYKFHFLPDLR\":0.999,\"ANSNDFVGNVDVDYKFHFLPDLR\":0.999,\"ANSNDFVGNVDVDYK\":0.999,\"FHFLPDLR\":0.999 }".to_string();
    let unipept_responses: String = "[{\"sequence\": \"AAAEAYFGK\", \"taxa\": [997888, 456163, 997892, 262724, 1235785, 226186, 762633, 338188, 702447, 483215, 371601, 818, 274, 28116, 798128, 300852, 869210, 657309]}, {\"sequence\": \"AADFLDGIK\", \"taxa\": [997888, 246787, 997892, 1796613, 1133319, 226186, 338188, 2447885, 483215, 371601, 763034, 2763675, 871324, 657309, 871325, 2650157, 818, 537012, 820, 2725562, 329854, 291644, 291645, 47678, 471870, 85831, 2838472, 742727, 1235785, 1235787, 2715212, 1841865, 28113, 1339346, 28116, 1339349, 411479, 693979, 674529, 1445607, 702447, 997873, 997874, 626931, 997884, 997886, 997887]}, {\"sequence\": \"AAGVSVQNVSGTFGTAPK\", \"taxa\": [679937, 470145, 246787, 997892, 1796613, 226186, 454154, 371601, 2608404, 679190, 310298, 2763675, 871324, 657309, 649761, 1115809, 2569763, 1236517, 1796646, 888743, 470565, 873513, 1703337, 1122985, 944557, 699437, 322095, 1703345, 818, 862515, 820, 537012, 700598, 823, 681398, 165179, 471870, 1263040, 947013, 2094150, 550983, 85831, 1235785, 1841864, 2838731, 2715212, 857291, 1235787, 997886, 1339346, 28116, 762968, 76123, 869213, 28127, 28129, 674529, 477666, 28132, 1220578, 2490854, 1177574, 762982, 2291815, 386414, 702447, 2992112, 1401073, 1401074, 2840690, 997874, 1602169, 2321403, 329854]}, {\"sequence\": \"AAIEEGIIPGGGVAYIR\", \"taxa\": [997888, 997892, 1796613, 226186, 454154, 338188, 483215, 371601, 2763675, 657309, 2650157, 818, 291644, 291645, 47678, 85831, 1235785, 2715212, 1339346, 28116, 762968, 674529, 762982, 702447, 997873, 1678841, 997884, 997886, 997887]}, {\"sequence\": \"AALENAQIIFDKK\", \"taxa\": [1547, 445974]}, {\"sequence\": \"AALQTSSFMSAASFQETTK\", \"taxa\": [997888, 246787, 997892, 1796613, 454154, 203275, 2447885, 387090, 310297, 310298, 310300, 1433126, 1796646, 397865, 2650157, 47678, 2094150, 2715212, 1099853, 2530390, 2530391, 762968, 2840674, 762982, 762984, 2840690, 329854, 470145, 449673, 2840716, 2840718, 763034, 272559, 2751153, 484018, 2590900, 2725562, 1349822, 1073351, 1841865, 2838732, 2838745, 693979, 674529, 547042, 1445607, 2838766, 626929, 862962, 626931, 2840839, 1133319, 338188, 2840856, 817, 818, 820, 291644, 291645, 471870, 742725, 742726, 742727, 85831, 1121097, 1235785, 1235787, 1121100, 411479, 1321819, 1235809, 665953, 667015, 226186, 483215, 1339280, 371601, 2763675, 871324, 871325, 657309, 46506, 1339311, 1339314, 1339315, 1339316, 537012, 681398, 1339327, 1263040, 2838468, 2838471, 2838472, 1339337, 2838474, 880074, 1263052, 28111, 28112, 28113, 1339346, 28116, 1339349, 477666, 28139, 295405, 702447, 997873, 997874, 997881, 997883, 997884, 997886, 997887]}]".to_string();
    let pep_psm_counts: String = "{\"GIVKDTTGEPVIGANVVVK\":9,\"DTTGEPVIGANVVVK\":4,\"GTTTGTITDFDGNFQLSAK\":14,\"KNDMTGSVMAIKPDELSK\":63,\"NDMTGSVMAIKPDELSK\":7,\"NDMTGSVMAIKPDELSKGITTNAQDMLSGK\":2,\"GITTNAQDMLSGK\":40,\"TTNAQDMLSGK\":1,\"IAGVSVISNDGTPGGGAQIR\":20,\"SVISNDGTPGGGAQIR\":1,\"ISNDGTPGGGAQIR\":1,\"SNDGTPGGGAQIR\":1,\"IRGGSSLNASNDPLIVIDGLAIDNEGIK\":12,\"GGSSLNASNDPLIVIDGLAIDNEGIK\":26,\"GMANGLSMVNPADIETLTVLK\":63,\"GMANGLSMVNPADIETLTVLKDASATAIYGSR\":5,\"DASATAIYGSR\":6,\"ASNGVIIITTK\":4,\"ASNGVIIITTKK\":3,\"SNGVIIITTKK\":1,\"NGQAPSVSYNGSVSFSK\":10,\"RYDVLSGDEYR\":13,\"YDVLSGDEYR\":4,\"AYANQLWGDK\":3,\"AYANQLWGDKLPADLGTANTDWQDQIFR\":9,\"LPADLGTANTDWQDQIFR\":17,\"TAVSTDHHVSINGGFK\":10,\"AVSTDHHVSINGGFK\":2,\"VSTDHHVSINGGFK\":1,\"STDHHVSINGGFK\":1,\"VSLGYTDDNGIVK\":4,\"FTASVNLAPSFFEDHLK\":44,\"ASVNLAPSFFEDHLK\":5,\"YADTGAAIGGALAIDPTRPVYSNEDPYQFTGGYWQNINSTTGFSNPDWK\":16,\"YTSNPNSPQNPLAALELKNDKANSNDFVGNVDVDYK\":19,\"YTSNPNSPQNPLAALELK\":11,\"YTSNPNSPQNPLAALELKNDK\":11,\"PNSPQNPLAALELK\":1,\"SPQNPLAALELK\":1,\"NDKANSNDFVGNVDVDYK\":9,\"NDKANSNDFVGNVDVDYKFHFLPDLR\":4,\"ANSNDFVGNVDVDYKFHFLPDLR\":5,\"ANSNDFVGNVDVDYK\":4,\"FHFLPDLR\":21,\"LHASIGGEYAEGTQTTIVSPYSFGNNYYGWNGDVTQYK\":21,\"SFGNNYYGWNGDVTQYK\":1,\"SLGANDFDIMVGGEEQHFHR\":101,\"NGFEEGQGWDSYTQEPHDAK\":15,\"GFEEGQGWDSYTQEPHDAK\":1,\"LREQTAYATR\":4,\"NTLVSYFGR\":6,\"LNYSLLNR\":2,\"IKEENFLKDVNVLSDLK\":22,\"IKEENFLKDVNVLSDLKLR\":1,\"IKEENFLK\":1,\"EENFLKDVNVLSDLK\":1,\"DVNVLSDLK\":4,\"LGWGITGQQNIGDDFAYLPLYVVNNEYAQYPFGDTYYSTSR\":2,\"LGWGITGQQNIGDDFAYLPLYVVNNEYAQYPFGDTYYSTSRPK\":2,\"AFNENLKWEK\":7,\"AFNENLK\":2,\"FNENLKWEK\":2,\"NENLKWEK\":2,\"TTTWNAGLDFGFLNGR\":11,\"ITGGIDGYFR\":13,\"ITGGIDGYFRK\":4,\"TGGIDGYFR\":1,\"KTDDLLNSVK\":3,\"TDDLLNSVK\":2,\"IPVGTNFNAQMTQNIGSLENYGMEFSINAKPIVTK\":10,\"IPVGTNFNAQMTQNIGSLENYGMEFSINAK\":3,\"ENYGMEFSINAKPIVTK\":2,\"DFTWDLSYNITWNHNEITK\":12,\"LTGGDDSDYYVEAGDKISR\":14,\"LTGGDDSDYYVEAGDK\":6,\"LTGGDDSDYYVEAGDKIS\":1,\"TGGDDSDYYVEAGDKISR\":1,\"GGDDSDYYVEAGDKISR\":1,\"VGYAANSFYVYQQVYDENGKPIENMFVDR\":25,\"NGNGTIDSGDKYIYK\":3,\"GNGTIDSGDKYIYK\":1,\"NGTIDSGDKYIYK\":2,\"KPAGDVLMGLTSK\":29,\"NFDFSFSLR\":17,\"ASLNNYVYYDFLSNK\":28,\"VYATVQNPFIISK\":11,\"VYATVQNPF\":2 }".to_string();

    perform_taxa_weighing(unipept_responses, pep_scores, pep_psm_counts, 10, "species".to_string());
}
