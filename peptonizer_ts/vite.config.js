import { defineConfig } from "vite";
import path from "path";

export default defineConfig({
    assetsInclude: ['**/*.py', '**/*.whl'],
    worker: {
        format: "es"
    },
    build: {
        lib: {
            entry: path.resolve(__dirname, 'src/index.ts'), // Change this to your library's entry point
            name: 'Peptonizer', // The global variable name for IIFE/UMD builds
            formats: ['es'], // Switch to ESM format,
            filename: "peptonizer.js"
        },
        rollupOptions: {
            // Externalize dependencies you don't want to bundle
            external: ['pyodide'], // Example, add others like 'react' if neededp
            output: {
                globals: {
                    pyodide: 'Pyodide', // Define the global name for external libraries
                }
            }
        },
    },
    optimizeDeps: { exclude: ["pyodide"] }
});
