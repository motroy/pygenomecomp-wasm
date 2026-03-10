// Service worker — caches Pyodide core files so repeat visits load instantly.
// The cache name is versioned to match Pyodide; bump it when upgrading Pyodide.
const CACHE_NAME = 'pyodide-v0.27.0';

const PYODIDE_URLS = [
  'https://cdn.jsdelivr.net/pyodide/v0.27.0/full/pyodide.js',
  'https://cdn.jsdelivr.net/pyodide/v0.27.0/full/pyodide.asm.wasm',
  'https://cdn.jsdelivr.net/pyodide/v0.27.0/full/python_stdlib.zip',
  'https://cdn.jsdelivr.net/pyodide/v0.27.0/full/pyodide-lock.json',
];

// Pre-cache Pyodide files when the SW is installed.
self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME).then(cache => cache.addAll(PYODIDE_URLS))
  );
  // Activate immediately — don't wait for existing tabs to close.
  self.skipWaiting();
});

// Remove old caches on activation.
self.addEventListener('activate', event => {
  event.waitUntil(
    caches.keys().then(keys =>
      Promise.all(keys.filter(k => k !== CACHE_NAME).map(k => caches.delete(k)))
    )
  );
  self.clients.claim();
});

// Cache-first for Pyodide URLs; network-only for everything else.
self.addEventListener('fetch', event => {
  const url = event.request.url;
  if (PYODIDE_URLS.some(u => url.startsWith(u.split('?')[0]))) {
    event.respondWith(
      caches.open(CACHE_NAME).then(cache =>
        cache.match(event.request).then(cached => {
          if (cached) return cached;
          return fetch(event.request).then(response => {
            cache.put(event.request, response.clone());
            return response;
          });
        })
      )
    );
  }
});
