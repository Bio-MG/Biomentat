#!/usr/bin/env node
// =============================================================================
// tools/ws_drive_client.mjs — P2.3 Option B : client WebSocket sans navigateur
// Maintient une session Shiny active pour piloter biomentat-drive sans Chrome.
// Usage : node tools/ws_drive_client.mjs --url=ws://127.0.0.1:6561/websocket/ [--timeout-s=30]
// =============================================================================

function parseArgs(argv) {
  const opts = { host: "127.0.0.1", port: 6561, url: null, timeoutS: 30 };
  for (const arg of argv) {
    if (arg.startsWith("--url=")) opts.url = arg.slice(6);
    else if (arg.startsWith("--host=")) opts.host = arg.slice(7);
    else if (arg.startsWith("--port=")) opts.port = Number(arg.slice(7)) || 6561;
    else if (arg.startsWith("--timeout-s=")) opts.timeoutS = Math.max(1, Number(arg.slice(12)) || 30);
  }
  if (!opts.url) opts.url = `ws://${opts.host}:${opts.port}/websocket/`;
  return opts;
}

const opts = parseArgs(process.argv.slice(2));
const ws = new WebSocket(opts.url);

const timer = setTimeout(() => {
  try { ws.close(1000, "timeout"); } catch (_) {}
  process.exit(0);
}, opts.timeoutS * 1000);

ws.addEventListener("open", () => {
  const initPayload = {
    method: "init",
    data: {
      ".clientdata_pixelratio": 1,
      ".clientdata_url_protocol": "http:",
      ".clientdata_url_hostname": opts.host,
      ".clientdata_url_port": String(opts.port),
      ".clientdata_url_pathname": "/",
      ".clientdata_url_search": "",
      ".clientdata_url_hash_initial": "",
      ".clientdata_singletons": ""
    }
  };
  ws.send(JSON.stringify(initPayload));
});

ws.addEventListener("message", (ev) => {
  let msg;
  try {
    msg = JSON.parse(String(ev.data));
  } catch (_) {
    return;
  }
  if (msg && msg.inputMessages && Array.isArray(msg.inputMessages)) {
    const data = {};
    for (const im of msg.inputMessages) {
      if (im && im.id && Object.prototype.hasOwnProperty.call(im.message || {}, "value")) {
        data[im.id] = im.message.value;
      }
    }
    if (Object.keys(data).length > 0) {
      ws.send(JSON.stringify({ method: "update", data }));
    }
  }
});

ws.addEventListener("error", (err) => {
  clearTimeout(timer);
  console.error("ws_drive_client error:", err?.message || err);
  process.exit(1);
});

ws.addEventListener("close", () => {
  clearTimeout(timer);
  process.exit(0);
});
