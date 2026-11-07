/* Jump to a physical page in the reader's local copy of the analysis note. */
(() => {
  "use strict";
  const views = new Set();
  let objectUrl = null;
  let request = 0;

  function refresh() {
    for (const view of views) {
      const page = Number(view.page());
      const source = objectUrl || (window.ATLAS_HAS_LOCAL_NOTE === true ? "./analysis-note.pdf" : null);
      if (source && Number.isInteger(page) && page > 0) {
        view.link.href = source + "#page=" + page;
        view.link.removeAttribute("aria-disabled");
        view.link.title = "Open page " + page + (objectUrl ? " of the PDF you selected" : " of the bundled analysis note");
      } else {
        view.link.removeAttribute("href");
        view.link.setAttribute("aria-disabled", "true");
        view.link.title = "Choose your local copy of the analysis note first";
      }
    }
  }

  function clearPdf() {
    if (objectUrl) URL.revokeObjectURL(objectUrl);
    objectUrl = null;
  }

  function attach({pick, link, page}) {
    link.target = "_blank";
    link.rel = "noopener";
    views.add({pick, link, page});
    pick.addEventListener("change", async () => {
      const activeRequest = ++request;
      clearPdf();
      refresh();
      pick.setCustomValidity("");
      const file = pick.files[0];
      if (!file) return;
      try {
        const header = await file.slice(0, 1024).text();
        if (activeRequest !== request) return;
        if (!header.includes("%PDF-")) throw new Error("Choose a PDF.");
        objectUrl = URL.createObjectURL(new Blob([file], {type: "application/pdf"}));
        refresh();
      } catch {
        if (activeRequest !== request) return;
        pick.setCustomValidity("Choose a readable PDF of the analysis note.");
        pick.reportValidity();
      }
    });
    refresh();
  }

  // Object URLs are released by the browser when this document unloads.
  window.notePage = {attach, refresh};
})();
