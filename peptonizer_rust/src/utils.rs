#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
use serde::Deserialize;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    pub fn log(s: &str);

    #[wasm_bindgen(js_namespace = console)]
    pub fn error(s: &str);
}

#[cfg(not(target_arch = "wasm32"))]
pub fn log(s: &str) {
    println!("{}", s);
}

#[cfg(not(target_arch = "wasm32"))]
pub fn error(s: &str) {
    eprintln!("{}", s);
}

#[derive(Debug, Deserialize, Clone)]
pub struct UnipeptJson {
    pub sequence: String,
    pub taxa: Vec<i32>,
}
