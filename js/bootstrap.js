// The home screen needs only the organism configuration. Atlas scripts run after a choice.
async function bootstrapAtlas() {
  const home = document.querySelector("#organism-home");
  if (!ORGANISM) {
    const cards = document.querySelector("#organism-cards");
    for (const organism of Object.values(ORGANISMS)) {
      const card = document.createElement("a");
      card.className = "organism-card";
      card.href = organismURL(organism.id);
      const name = document.createElement("h2");
      name.textContent = organism.name;
      const scientificName = document.createElement("p");
      scientificName.className = "home-scientific-name";
      scientificName.textContent = organism.scientificName;
      const action = document.createElement("span");
      action.className = "home-open-label";
      action.textContent = "Open atlas →";
      card.append(name, scientificName, action);
      cards.append(card);
    }
    home.hidden = false;
    return;
  }

  document.querySelector("#atlas-app").hidden = false;
  initializeOrganismUI();

  // Preserve the existing classic-script order and shared globals; app.js starts initialization.
  const scripts = [
    "protein-annotations", "analysis", "control-qc", "go-analysis", "compact-data",
    "ligand-viewer", "pocket-electrostatics", "pose-electrostatics", "pocket-cloud",
    "protein-viewer", "app",
  ];
  try {
    for (const name of scripts) {
      if (organismLeaving) return;
      await new Promise((resolve, reject) => {
        const script = document.createElement("script");
        script.src = `js/${name}.js?v=yeast-organism-1`;
        script.async = false;
        script.onload = resolve;
        script.onerror = () => reject(new Error(`Could not load js/${name}.js. Refresh the page to retry.`));
        document.body.append(script);
      });
    }
  } catch (error) {
    if (organismLeaving) return;
    const title = document.createElement("h2");
    title.textContent = "Atlas could not start";
    const message = document.createElement("p");
    message.textContent = error.message;
    document.querySelector("#loading-screen").replaceChildren(title, message);
  }
}

bootstrapAtlas();
