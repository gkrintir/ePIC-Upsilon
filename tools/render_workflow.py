#!/usr/bin/env python3
"""Render the README workflow graphic using only the standard library."""
from pathlib import Path
from html import escape

ROOT=Path(__file__).resolve().parents[1]
W,H=1080,910
nodes=[
 ('inputs',55,25,440,76,'1  Reconstructed ePIC events','Download / explicit input list','external'),
 ('existing',585,25,440,76,'Existing merged UpsilonTree','Start here when tuples are already available','external'),
 ('produce',55,145,440,76,'2  ROOT producer + Condor','Decay electrons, scattered electron, recoil proton','production'),
 ('merge',55,265,440,76,'3  Merge job outputs','Retain generated-event counts and provenance','production'),
 ('normalize',320,385,440,76,'4  Normalize and audit','State rate, luminosity, branches and photon flux','analysis'),
 ('full',55,505,440,86,'5  Full state analysis','Mass, efficiency, cross sections and resolutions','analysis'),
 ('opt',585,505,440,86,'6  Reoptimization','Response, binning, inference and stress tests','analysis'),
 ('review',320,635,440,76,'7  Review validation status','Accepted conditional results / diagnostics / unavailable','review'),
 ('plot',320,745,440,76,'8  SWAN notebooks and exports','One profile selects runs, category and luminosity','plot'),
]
colors={'external':('#f1f5f9','#94a3b8'),'production':('#fef3c7','#d97706'),
        'analysis':('#e0f2fe','#0284c7'),'review':('#fce7f3','#be185d'),'plot':('#dcfce7','#16a34a')}
paths=['M275 101 V145','M275 221 V265','M275 341 V365 H540 V385',
       'M805 101 V365 H540 V385','M540 461 V483 H275 V505','M540 461 V483 H805 V505',
       'M275 591 V615 H540 V635','M805 591 V615 H540 V635','M540 711 V745']
parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title desc">',
 '<title id="title">ePIC Upsilon workflow from event files to reviewed figures</title>',
 '<desc id="desc">New production and existing merged inputs join at normalization. Full analysis and reoptimization run separately before validation and plotting.</desc>',
 '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#475569"/></marker></defs>',
 '<rect width="1080" height="910" rx="14" fill="#ffffff"/>']
for d in paths:parts.append(f'<path d="{d}" fill="none" stroke="#475569" stroke-width="2" marker-end="url(#arrow)"/>')
for name,x,y,w,h,title,sub,kind in nodes:
 fill,stroke=colors[kind]
 parts.extend([f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>',
 f'<text x="{x+20}" y="{y+30}" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#0f172a">{escape(title)}</text>',
 f'<text x="{x+20}" y="{y+55}" font-family="Arial,sans-serif" font-size="14" fill="#334155">{escape(sub)}</text>'])
parts.extend(['<text x="540" y="858" text-anchor="middle" font-family="Arial,sans-serif" font-size="17" fill="#0f172a">Separate inclusive / proton-tagged workspaces · EOS data stay external</text>',
 '<text x="540" y="886" text-anchor="middle" font-family="Arial,sans-serif" font-size="14" fill="#475569">Simulation cards are pinned references. Production requires the documented source and environment checks.</text>','</svg>'])
(ROOT/'docs/workflow.svg').write_text('\n'.join(parts)+'\n')
