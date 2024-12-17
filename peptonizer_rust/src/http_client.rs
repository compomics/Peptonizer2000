use serde::{Serialize, Deserialize };

#[cfg(target_arch = "wasm32")]
use web_sys::{XmlHttpRequest};

pub trait HttpClient {
    fn perform_post_request(&self, url: String, batch: Vec<i32>) -> Result<String, String>;
}

#[cfg(target_arch = "wasm32")]
pub struct WasmHttpClient;

#[cfg(target_arch = "wasm32")]
impl HttpClient for WasmHttpClient {
    fn perform_post_request(&self, url: String, batch: Vec<i32>) -> Result<String, String> {
        let payload = HTTPPostPayload {
            input: batch,
            extra: true
        };
        let payload_json = serde_json::to_string(&payload).unwrap();

        // Create a new XMLHttpRequest object
        let xhr = XmlHttpRequest::new().map_err(|e| format!("Failed to create XMLHttpRequest: {:?}", e))?;
        
        // Open the request (synchronous mode by setting `async` to false)
        xhr.open_with_async("POST", &url, false)
            .map_err(|e| format!("Failed to open request: {:?}", e))?;
        
        // Set the request header for JSON
        xhr.set_request_header("Content-Type", "application/json")
            .map_err(|e| format!("Failed to set request header: {:?}", e))?;
        
        // Send the request with the body
        xhr.send_with_opt_str(Some(&payload_json))
            .map_err(|e| format!("Failed to send request: {:?}", e))?;
        
        let status = xhr.status().map_err(|_e| format!("Failed to extract status from response"))?;
        if status == 200 {
            let response = xhr.response_text()
            .expect("Expected json in response")
            .ok_or(format!("Failed to extract text from response"))?;
        
            return Ok(format!("{}", response));
        }

        Err(format!("Status code {}", status))
    }
}

pub fn create_http_client() -> impl HttpClient {
    #[cfg(target_arch = "wasm32")]
    WasmHttpClient
}


#[derive(Serialize, Deserialize, Debug)]
struct HTTPPostPayload {
    input: Vec<i32>,
    extra: bool
}