[ {"selector": "node", "css": {
    "shape": "ellipse",
    "text-valign":"center",
    "text-halign":"center",
    "content": "data(id)",
    "background-color": "red",
    "border-color": "black","border-width":"1px",
    "width":  "40px",
    "height": "40px",
    "font-size":"14px"}},

 {"selector":"node:selected", "css": {
     "text-valign":"center",
     "text-halign":"center",
     "border-color": "black",
     "overlay-opacity": 0.2,
     "overlay-color": "gray"
     }},
     

  {"selector": "node[[degree = 1]]", "css": {
       "background-color": "pink"
  }},

  {"selector": "node[[degree >= 2]][[degree < 5]]", "css": {
      "width": "50px",
      "height": "50px"
  }},

  {"selector": "node[[degree >= 5]][[degree < 8]]", "css": {
      "width": "60px",
      "height": "60px"
  }},

  {"selector": "node[[degree >= 8]][[degree < 12]]", "css": {
      "width": "70px",
      "height": "70px"
  }},

  {"selector": "node[[degree >= 12]][[degree < 17]]", "css": {
      "width": "80px",
      "height": "80px"
  }},

  {"selector": "node[[degree >= 17]]", "css": {
      "width": "100px",
      "height": "100px"
  }},

 {"selector": "edge", "css": {
     "opacity": 0.5,
     "curve-style": "bezier"
     }},

 {"selector":"edge:selected", "css": {
     "overlay-opacity": 0.2,
     "overlay-color": "red"
     }}
]
