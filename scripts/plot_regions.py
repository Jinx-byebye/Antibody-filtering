#!/usr/bin/env python3
"""Draw the matched regional-confidence study from evaluated numeric results."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter
import numpy as np

METRICS = ("nac", "n3ac", "nac_cdr3", "antibody", "cdr", "cdr3", "interface")
LABELS = ("NAC", r"N$^3$AC", "CDR3+NAC", "Whole antibody", "CDR", "CDR3 only", "Interface")
PREDICTORS = (("af3", "AlphaFold 3"), ("protenix", "Protenix-v2"), ("esmfold2", "ESMFold2"))
COLORS = {"ngs": "#266B8D", "igdesign": "#C56B2C"}
INK = "#233640"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    p.add_argument("--output", type=Path, required=True, help="Output stem")
    args = p.parse_args()
    records = list(csv.DictReader((args.results / "summary.csv").open()))
    predictors = PREDICTORS + ((("opendde", "OpenDDE"),) if any(r["predictor"] == "opendde" for r in records) else ())
    assert len(records) == len(predictors)*42
    data = {(r["predictor"],r["dataset"],r["metric"],r["endpoint"]):
            tuple(float(r[k]) for k in ("mean","low","high")) for r in records}
    assert len(data) == len(predictors)*42
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10.5,
        "mathtext.fontset":"dejavusans", "pdf.fonttype":42, "ps.fonttype":42,
        "svg.fonttype":"none", "axes.labelcolor":INK, "text.color":INK,
        "xtick.color":INK, "ytick.color":INK, "axes.unicode_minus":False})
    fig = plt.figure(figsize=(7.2,9.9 if len(predictors)==4 else 9.3),facecolor="white")
    fig.text(.025,.978,"A",weight="bold",size=13)
    fig.text(.065,.978,"Which confidence is averaged?",weight="bold",size=12)
    fig.text(.065,.955,"Full antibody + antigen are co-folded for every score",size=10.5,color="#536873")
    mask_ax = fig.add_axes([.025,.778,.95,.163])
    mask_ax.set_xlim(0,1); mask_ax.set_ylim(-.6,7.7); mask_ax.axis("off")
    regions=("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")
    include=((1,0,1,0,1,0,1),(1,1,1,1,1,0,1),(1,0,1,0,1,1,1),
             (1,1,1,1,1,1,1),(0,1,0,1,0,1,0),(0,0,0,0,0,1,0))
    x0,step=.245,.073
    for j,name in enumerate(regions):
        mask_ax.text(x0+step*j+.031,7.3,name,ha="center",va="center",size=9.9,
                     color=INK if j%2==0 else "#956334")
    for i,label in enumerate(LABELS):
        y=6.4-i
        mask_ax.text(.01,y,label,va="center",weight="bold" if i==0 else "normal",size=10.5)
        if i < 6:
            for j,on in enumerate(include[i]):
                mask_ax.add_patch(Rectangle((x0+step*j,y-.27),.062,.54,
                    facecolor="#758994" if on else "white",edgecolor="#C7D0D5",linewidth=.65))
        else:
            mask_ax.text(x0,y,"Contact residues on both partners",va="center",size=10.3)
    mask_ax.add_patch(Rectangle((.80,5.95),.035,.5,facecolor="#758994",linewidth=0))
    mask_ax.text(.85,6.2,"Included",va="center",size=10)
    mask_ax.add_patch(Rectangle((.80,4.95),.035,.5,facecolor="white",edgecolor="#C7D0D5",linewidth=.65))
    mask_ax.text(.85,5.2,"Excluded",va="center",size=10)
    mask_ax.text(.80,3.6,"Same rule for\nheavy and\nlight chains",va="top",size=10,linespacing=1.4,color="#536873")
    mask_ax.text(.80,.0,"8 Å contacts",va="center",size=10,color="#536873")
    fig.add_artist(Line2D([.025,.975],[.766,.766],transform=fig.transFigure,color="#D6DEE2",lw=.8))
    fig.text(.025,.742,"B",weight="bold",size=13)
    fig.text(.065,.742,"From binder separability to screening gains",weight="bold",size=12)
    handles=[Line2D([],[],color=COLORS["ngs"],marker="o",ls="none",ms=5,label="Display-based enrichment"),
             Line2D([],[],color=COLORS["igdesign"],marker="D",ls="none",ms=4.5,label="IgDesign SPR")]
    fig.legend(handles=handles,loc="center left",bbox_to_anchor=(.055,.714),ncol=2,
               frameon=False,columnspacing=1.7,handletextpad=.4,fontsize=10.5)
    xs=(.225,.491,.757); width=.213
    endpoints=("auc","recall_q_20","utility")
    headers=("Separability","Fixed budget","Best threshold")
    subheaders=("AUC",r"Recall at $Q=20\%$",r"Utility $U^\star$")
    for x,h,s in zip(xs,headers,subheaders):
        fig.text(x+width/2,.682,h,ha="center",weight="bold",size=10.7)
        fig.text(x+width/2,.663,s,ha="center",size=10)
    bounds=[]
    for endpoint in endpoints:
        values=[v for k,v in data.items() if k[-1]==endpoint]
        lows=[v[1] for v in values]; highs=[v[2] for v in values]
        if endpoint == "auc":
            bounds.append((max(0,np.floor((min(lows)-.025)*10)/10),min(1,np.ceil((max(highs)+.025)*10)/10)))
        else:
            bounds.append((0,min(1,np.ceil((max(highs)+.025)*10)/10)))
    for i,(pred,title) in enumerate(predictors):
        bottom=.494-i*.153 if len(predictors)==4 else .452-i*.215
        height=.120 if len(predictors)==4 else .173
        fig.text(.025,bottom+height+.013,title,weight="bold",size=11.3)
        for j,endpoint in enumerate(endpoints):
            ax=fig.add_axes([xs[j],bottom,width,height])
            ax.axhspan(-.48,.48,color="#EDF3F5",zorder=0)
            for d,(ds,marker) in enumerate((("ngs","o"),("igdesign","D"))):
                y=np.arange(7)+(d-.5)*.29
                values=np.array([data[pred,ds,m,endpoint] for m in METRICS])
                mean,lo,hi=values.T
                ax.hlines(y,lo,hi,color=COLORS[ds],lw=1.0,alpha=.7,zorder=2)
                ax.scatter(mean,y,c=COLORS[ds],marker=marker,s=22 if ds=="ngs" else 19,
                           edgecolor="white",linewidth=.4,zorder=3)
            if j<2:
                ax.axvline(.5 if j==0 else .2,color="#8C9CA5",lw=.8,ls=(0,(3,3)),zorder=1)
            ax.set_xlim(*bounds[j]); ax.set_ylim(6.65,-.65)
            ax.set_yticks(range(7))
            ax.set_yticklabels(LABELS if j==0 else [],fontsize=10.5)
            if j==0:
                ax.get_yticklabels()[0].set_fontweight("bold")
            ax.tick_params(axis="y",length=0,pad=8)
            ax.tick_params(axis="x",length=3,width=.6,pad=3,labelsize=10)
            ticks = ((.5,.7,.9),(.2,.5,.8),(0,.2,.4,.6))[j]
            ax.set_xticks([t for t in ticks if bounds[j][0] <= t <= bounds[j][1]])
            if j==1:
                ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
            ax.grid(axis="x",color="#E7EBEE",lw=.5,zorder=0)
            for side in ("top","right","left"):
                ax.spines[side].set_visible(False)
            ax.spines["bottom"].set_color("#B6C2C9")
            ax.spines["bottom"].set_linewidth(.6)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","svg","png"):
        fig.savefig(args.output.with_suffix("."+ext),dpi=300,facecolor="white")
    (args.results/"figure_data.json").write_text(json.dumps(records,indent=2))
    plt.close(fig)


if __name__=="__main__":
    main()
