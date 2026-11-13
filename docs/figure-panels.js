/* These view boxes select the original plot artwork. No values are read from pixels. */
(() => {
  const dimensions = window.PAPER?.image_dimensions || {mass: [1600,1424], w:[1477,1600], differential:[1486,1600]};
  const boxes = {
    mass: {mass: ["(a) Three-state mass fit",[0,0,1,.53]]},
    efficiency: {pt_y:["(b) Transverse momentum and rapidity",[0,.514,.506,.486]],xv_w:["(c) xᵥ and W",[.50,.514,.50,.486]]},
    w: {"1S":["(a) Υ(1S)",[0,0,1,.56]],"2S":["(b) Υ(2S)",[0,.545,.525,.455]],"3S":["(c) Υ(3S)",[.505,.545,.495,.455]]},
    differential: {pt2:["(a) Dielectron pT² and slope",[0,0,1,.56]],y:["(b) Laboratory rapidity",[0,.55,1,.45]]}
  };
  const ns="http://www.w3.org/2000/svg";
  window.showFigurePanel = function({view,image,src,selection,onChange}) {
    const select=document.getElementById("plot-panel");
    const full=document.getElementById("figure-image");
    const area=document.getElementById("panel-zoom");
    select.replaceChildren(new Option("Full paper figure","all"));
    for(const [key,[label]] of Object.entries(boxes[view])) select.add(new Option(label,key));
    select.value=selection in boxes[view] ? selection : "all";
    function draw() {
      area.replaceChildren();
      const key=select.value;
      if(key==="all") {area.style.display="none";full.style.display="block";return;}
      const [width,height]=dimensions[image];
      const [label,[x,y,w,h]]=boxes[view][key];
      area.setAttribute("viewBox",[x*width,y*height,w*width,h*height].join(" "));
      area.setAttribute("aria-label",label);
      area.style.maxWidth=(440*w*width/(h*height))+"px";
      area.style.height="auto";
      area.style.margin="0 auto";
      const art=document.createElementNS(ns,"image");
      art.setAttribute("width",width);art.setAttribute("height",height);art.setAttribute("href",src);
      area.append(art);area.style.display="block";full.style.display="none";
    }
    select.onchange=()=>{draw();onChange(select.value);};
    draw();
  };
})();
