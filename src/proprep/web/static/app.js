/* ProPrep Web Shell — terminal pane.
 *
 * Connects an xterm.js terminal to the /ws/term websocket. Server frames
 * are bytes (PTY output) or text JSON control messages; client frames are
 * bytes (stdin) or text JSON for resize.
 */

(function () {
  "use strict";

  const statusEl = document.getElementById("status");
  const termHost = document.getElementById("term");

  function setStatus(text, kind) {
    statusEl.textContent = text;
    statusEl.classList.remove("ok", "err");
    if (kind) statusEl.classList.add(kind);
  }

  // Full 16-color ANSI palettes (not just bg/fg/cursor) so Rich's *named*
  // colors render the same here as in a native terminal instead of falling
  // back to xterm's built-in defaults. The light palette uses darker hues so
  // they stay legible on white. Truecolor hex from Rich (e.g. the protonation
  // determinant green/orange) bypasses these and renders exactly.
  const THEMES = {
    dark: {
      background: "#000000", foreground: "#e6edf3", cursor: "#58a6ff",
      selectionBackground: "#3b5070",
      black: "#484f58", red: "#ff7b72", green: "#3fb950", yellow: "#d29922",
      blue: "#58a6ff", magenta: "#bc8cff", cyan: "#39c5cf", white: "#b1bac4",
      brightBlack: "#6e7681", brightRed: "#ffa198", brightGreen: "#56d364",
      brightYellow: "#e3b341", brightBlue: "#79c0ff", brightMagenta: "#d2a8ff",
      brightCyan: "#56d4dd", brightWhite: "#f0f6fc",
    },
    light: {
      background: "#ffffff", foreground: "#24292f", cursor: "#0969da",
      selectionBackground: "#b6d7ff",
      black: "#24292f", red: "#cf222e", green: "#116329", yellow: "#7d4e00",
      blue: "#0969da", magenta: "#8250df", cyan: "#1b7c83", white: "#6e7781",
      brightBlack: "#57606a", brightRed: "#a40e26", brightGreen: "#1a7f37",
      brightYellow: "#633c01", brightBlue: "#218bff", brightMagenta: "#a475f9",
      brightCyan: "#3192aa", brightWhite: "#8c959f",
    },
  };

  // Terminal font size. Kept in localStorage so a size chosen for, say, a
  // manuscript screenshot survives reloads; proprep-web --font-size overrides
  // it at launch (see the /shell-theme fetch below).
  const FONT_KEY = "proprep.web.fontSize";
  const FONT_DEFAULT = 13, FONT_MIN = 8, FONT_MAX = 32;

  function clampFont(px) {
    if (!Number.isFinite(px)) return FONT_DEFAULT;
    return Math.min(FONT_MAX, Math.max(FONT_MIN, Math.round(px)));
  }

  function storedFontSize() {
    try {
      const raw = localStorage.getItem(FONT_KEY);
      if (raw === null) return FONT_DEFAULT;
      return clampFont(parseInt(raw, 10));
    } catch (_) {
      return FONT_DEFAULT;  // storage blocked (private mode / sandboxed frame)
    }
  }

  const term = new Terminal({
    cursorBlink: true,
    fontFamily: 'Menlo, Consolas, "DejaVu Sans Mono", monospace',
    fontSize: storedFontSize(),
    theme: THEMES.dark,
    // Lines retained above the viewport. ProPrep sessions are long, so this is
    // generous; xterm keeps scrollback in browser memory only, so it is still
    // ephemeral (a reload clears it). For a durable copy of the whole session,
    // the server also tees output to ~/.proprep/web_sessions/ — see
    // pty_session.py. Bumped from 5000 (2026-07).
    scrollback: 50000,
    convertEol: false,
    allowProposedApi: true,
  });
  const fitAddon = new FitAddon.FitAddon();
  term.loadAddon(fitAddon);
  term.loadAddon(new WebLinksAddon.WebLinksAddon());
  term.open(termHost);

  // Apply the launch-time theme (proprep-web --theme light|dark
  // [--high-contrast]). Until this resolves we render the dark default; on
  // failure we keep it. high-contrast uses xterm's minimumContrastRatio=7,
  // which auto-bumps any low-contrast text (incl. truecolor) against the bg.
  fetch("/shell-theme")
    .then((r) => r.json())
    .then((cfg) => {
      const name = cfg.theme === "light" ? "light" : "dark";
      term.options.theme = THEMES[name];
      term.options.minimumContrastRatio = cfg.highContrast ? 7 : 1;
      document.body.dataset.theme = name;
      // Only when --font-size was actually passed; otherwise the remembered
      // in-page size stands.
      if (typeof cfg.fontSize === "number") setFontSize(cfg.fontSize);
    })
    .catch(() => { /* keep the dark default */ });

  function fit() {
    try { fitAddon.fit(); } catch (_) { /* no-op pre-attach */ }
  }
  fit();

  // ---- websocket ---------------------------------------------------------

  const wsUrl = (location.protocol === "https:" ? "wss:" : "ws:") + "//" + location.host + "/ws/term";
  let ws = null;
  let reconnectTimer = null;
  let lastSentSize = { cols: 0, rows: 0 };

  function sendResize() {
    if (!ws || ws.readyState !== WebSocket.OPEN) return;
    const cols = term.cols, rows = term.rows;
    if (cols === lastSentSize.cols && rows === lastSentSize.rows) return;
    lastSentSize = { cols, rows };
    ws.send(JSON.stringify({ type: "resize", cols, rows }));
  }

  function connect() {
    setStatus("connecting…");
    ws = new WebSocket(wsUrl);
    ws.binaryType = "arraybuffer";

    ws.onopen = () => {
      setStatus("connected", "ok");
      // Push current size after the server has started the PTY.
      lastSentSize = { cols: 0, rows: 0 };
      sendResize();
      term.focus();
    };

    ws.onmessage = (ev) => {
      if (typeof ev.data === "string") {
        let msg;
        try { msg = JSON.parse(ev.data); } catch (_) { return; }
        if (msg.type === "exit") {
          term.write(`\r\n\x1b[33m[proprep exited with code ${msg.returncode}]\x1b[0m\r\n`);
          setStatus("exited (" + msg.returncode + ")", "err");
        } else if (msg.type === "error") {
          term.write(`\r\n\x1b[31m[shell error: ${msg.message}]\x1b[0m\r\n`);
          setStatus("error", "err");
        }
        return;
      }
      // Binary frame from PTY.
      const bytes = ev.data instanceof ArrayBuffer ? new Uint8Array(ev.data) : ev.data;
      term.write(bytes);
    };

    ws.onclose = () => {
      setStatus("disconnected", "err");
      // No automatic reconnect in Phase A — a closed shell means the child
      // exited; reloading the page is the explicit gesture to start fresh.
    };

    ws.onerror = () => {
      setStatus("ws error", "err");
    };
  }

  // ---- input + resize wiring --------------------------------------------

  term.onData((data) => {
    if (!ws || ws.readyState !== WebSocket.OPEN) return;
    // xterm.js gives us a JS string; the PTY child wants bytes. Encode UTF-8.
    ws.send(new TextEncoder().encode(data));
  });

  // Refit on window resize and after fonts settle.
  let resizeRaf = 0;
  function scheduleFit() {
    if (resizeRaf) cancelAnimationFrame(resizeRaf);
    resizeRaf = requestAnimationFrame(() => {
      resizeRaf = 0;
      fit();
      sendResize();
    });
  }
  window.addEventListener("resize", scheduleFit);
  // ResizeObserver picks up splitter drags too.
  if (window.ResizeObserver) {
    new ResizeObserver(scheduleFit).observe(termHost);
  }
  // After fonts/layout settle, refit once more.
  window.addEventListener("load", () => setTimeout(scheduleFit, 50));

  // ---- font size --------------------------------------------------------

  const fontReadout = document.getElementById("font-readout");
  const fontDec = document.getElementById("font-dec");
  const fontInc = document.getElementById("font-inc");

  function setFontSize(px, persist) {
    px = clampFont(px);
    // xterm 5 applies option writes live, so no reopen is needed. Bigger glyphs
    // mean fewer columns, so refit and tell the PTY — otherwise Rich keeps
    // wrapping tables to the old width.
    term.options.fontSize = px;
    fontReadout.textContent = px + "px";
    fontDec.disabled = px <= FONT_MIN;
    fontInc.disabled = px >= FONT_MAX;
    if (persist !== false) {
      try { localStorage.setItem(FONT_KEY, String(px)); } catch (_) { /* storage blocked */ }
    }
    scheduleFit();
  }

  fontDec.addEventListener("click", () => setFontSize(term.options.fontSize - 1));
  fontInc.addEventListener("click", () => setFontSize(term.options.fontSize + 1));

  // Ctrl/Cmd +/-/0, the shortcut every terminal emulator uses. We take it over
  // from the browser's page zoom on purpose: zoom would scale the whole UI and
  // leave xterm re-rasterizing at a fractional device pixel ratio, which is
  // exactly the blur to avoid in a figure.
  function fontKeyDelta(e) {
    if (!(e.ctrlKey || e.metaKey) || e.altKey) return null;
    if (e.key === "+" || e.key === "=") return 1;
    if (e.key === "-" || e.key === "_") return -1;
    if (e.key === "0") return 0;
    return null;
  }

  document.addEventListener("keydown", (e) => {
    const delta = fontKeyDelta(e);
    if (delta === null) return;
    e.preventDefault();
    setFontSize(delta === 0 ? FONT_DEFAULT : term.options.fontSize + delta);
  });

  // Returning false stops xterm from consuming the chord and sending it to the
  // PTY; the event still bubbles to the document listener above.
  term.attachCustomKeyEventHandler((e) => {
    if (e.type === "keydown" && fontKeyDelta(e) !== null) return false;
    return true;
  });

  // Sync the readout/buttons with the size the terminal was constructed with.
  setFontSize(term.options.fontSize, false);

  // ---- splitter ---------------------------------------------------------

  const layout = document.getElementById("layout");
  const splitter = document.getElementById("splitter");
  const FRAC_KEY = "proprep.web.termFrac";
  const MIN_FRAC = 0.15, MAX_FRAC = 0.85;

  function setFrac(frac) {
    frac = Math.min(MAX_FRAC, Math.max(MIN_FRAC, frac));
    layout.style.setProperty("--term-frac", String(frac));
    return frac;
  }

  // Restore saved width.
  const saved = parseFloat(localStorage.getItem(FRAC_KEY) || "");
  if (!Number.isNaN(saved)) setFrac(saved);

  function startDrag(startEvt) {
    startEvt.preventDefault();
    const rect = layout.getBoundingClientRect();
    const splitterW = splitter.getBoundingClientRect().width;
    layout.classList.add("dragging");

    function onMove(ev) {
      const x = (ev.touches ? ev.touches[0].clientX : ev.clientX) - rect.left;
      // Avoid divide-by-zero on degenerate widths.
      const usable = Math.max(1, rect.width - splitterW);
      setFrac(x / usable);
    }
    function onUp() {
      layout.classList.remove("dragging");
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
      window.removeEventListener("touchmove", onMove);
      window.removeEventListener("touchend", onUp);
      const cur = parseFloat(getComputedStyle(layout).getPropertyValue("--term-frac"));
      if (!Number.isNaN(cur)) localStorage.setItem(FRAC_KEY, String(cur));
    }
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
    window.addEventListener("touchmove", onMove, { passive: false });
    window.addEventListener("touchend", onUp);
  }
  splitter.addEventListener("mousedown", startDrag);
  splitter.addEventListener("touchstart", startDrag, { passive: false });
  // Keyboard nudges (focus the bar, then arrow keys).
  splitter.addEventListener("keydown", (ev) => {
    const cur = parseFloat(getComputedStyle(layout).getPropertyValue("--term-frac")) || 0.5;
    const step = ev.shiftKey ? 0.05 : 0.02;
    if (ev.key === "ArrowLeft")       { ev.preventDefault(); localStorage.setItem(FRAC_KEY, String(setFrac(cur - step))); }
    else if (ev.key === "ArrowRight") { ev.preventDefault(); localStorage.setItem(FRAC_KEY, String(setFrac(cur + step))); }
  });

  connect();

  // ---- viewer iframe / detached window via control websocket ------------

  const viewerFrame = document.getElementById("viewer-frame");
  const viewerPlaceholder = document.getElementById("viewer-placeholder");
  const detachBtn = document.getElementById("detach-btn");
  let lastAnnouncedPort = null;
  let detachedWindow = null;
  let detachedWatchdog = 0;

  function viewerUrl(port) {
    return "/viewer?_v=" + encodeURIComponent(port) + "&t=" + Date.now();
  }

  function showViewer(port) {
    // Route the new viewer URL to whichever target is currently active —
    // the docked iframe or the detached popup. Cache-bust on each
    // (re)launch so a server restart on the same port forces a fresh
    // load rather than reusing stale state.
    lastAnnouncedPort = port;
    const url = viewerUrl(port);
    if (detachedWindow && !detachedWindow.closed) {
      try { detachedWindow.location.replace(url); } catch (_) { /* opener gone */ }
      return;
    }
    viewerFrame.src = url;
    viewerFrame.hidden = false;
    if (viewerPlaceholder) viewerPlaceholder.hidden = true;
  }

  // ---- detach / redock --------------------------------------------------

  function applyDetachedLayout(detached) {
    layout.classList.toggle("detached", detached);
    detachBtn.textContent = detached ? "Re-dock viewer" : "Detach viewer";
    // Refit xterm to the new column count.
    scheduleFit();
  }

  function detach() {
    // Open (or re-focus) the popup. ``noopener`` would let the popup
    // outlive us, but we want to ``.close()`` it on re-dock — keep the
    // opener relationship.
    const url = lastAnnouncedPort ? viewerUrl(lastAnnouncedPort) : "/viewer";
    const features = "popup=yes,width=900,height=720";
    detachedWindow = window.open(url, "proprep-viewer", features);
    if (!detachedWindow) {
      // Pop-up blocker; surface a status hint.
      setStatus("popup blocked", "err");
      detachedWindow = null;
      return;
    }
    // Stop the docked iframe from polling while detached so we don't
    // run two pollers against the same backend.
    viewerFrame.src = "about:blank";
    viewerFrame.hidden = true;
    if (viewerPlaceholder) viewerPlaceholder.hidden = false;
    applyDetachedLayout(true);

    // Watch for the user closing the popup manually — revert state if so.
    clearInterval(detachedWatchdog);
    detachedWatchdog = setInterval(() => {
      if (!detachedWindow || detachedWindow.closed) {
        clearInterval(detachedWatchdog);
        detachedWatchdog = 0;
        detachedWindow = null;
        // Restore the docked iframe with the latest known viewer state.
        applyDetachedLayout(false);
        if (lastAnnouncedPort != null) showViewer(lastAnnouncedPort);
      }
    }, 500);
  }

  function redock() {
    if (detachedWindow && !detachedWindow.closed) {
      try { detachedWindow.close(); } catch (_) { /* ok */ }
    }
    clearInterval(detachedWatchdog);
    detachedWatchdog = 0;
    detachedWindow = null;
    applyDetachedLayout(false);
    if (lastAnnouncedPort != null) showViewer(lastAnnouncedPort);
  }

  detachBtn.addEventListener("click", () => {
    const isDetached = layout.classList.contains("detached");
    if (isDetached) redock(); else detach();
  });

  // If the parent page itself is closed, take the popup down with it so
  // a stale standalone window doesn't outlive the shell.
  window.addEventListener("pagehide", () => {
    if (detachedWindow && !detachedWindow.closed) {
      try { detachedWindow.close(); } catch (_) { /* ok */ }
    }
  });

  const controlWsUrl = (location.protocol === "https:" ? "wss:" : "ws:") + "//" + location.host + "/ws/control";
  let controlWs = null;
  let controlReconnect = null;

  function connectControl() {
    controlWs = new WebSocket(controlWsUrl);
    controlWs.onopen = () => { /* server replays cached state */ };
    controlWs.onmessage = (ev) => {
      let msg;
      try { msg = JSON.parse(ev.data); } catch (_) { return; }
      if (msg.type === "viewer_announce" && Number.isFinite(msg.port)) {
        // Always re-source — even if the port matches a previous
        // announcement, this catches "viewer was relaunched on same port".
        showViewer(msg.port);
      }
    };
    controlWs.onclose = () => {
      // The terminal websocket carries the auto-shutdown signal, not this
      // one. Reconnect quietly so a brief network blip doesn't leave the
      // viewer pane stuck.
      clearTimeout(controlReconnect);
      controlReconnect = setTimeout(connectControl, 1000);
    };
    controlWs.onerror = () => { /* close handler will retry */ };
  }
  connectControl();
})();
