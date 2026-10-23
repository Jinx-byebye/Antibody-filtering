#!/usr/bin/env python3
"""Render the reference tables as readable README figures."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from abfilter.tables import NAMES


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,default=Path('reference'))
    p.add_argument('--output',type=Path,default=Path('assets'))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    for table,stem,title in [('table2','across-antigen-results','Across-antigen filtering'),
                             ('table3','within-antigen-results','Within-antigen filtering')]:
        rows=list(csv.DictReader((a.reference/f'{table}.csv').open()))
        predictors=list(dict.fromkeys(r['predictor'] for r in rows))
        columns=list(rows[0])[2:]
        labels=['Metric','Recall\nFPR 10%','Recall\nFPR 20%','Recall\nFPR 30%',
                'Recall\nQ = 5%','Recall\nQ = 20%','Utility','Recall\nat U*','Precision\nat U*','Average\nrank']
        fig=plt.figure(figsize=(15,16.6),facecolor='white')
        fig.text(.02,.983,title,size=20,weight='bold',color='#263D48')
        fig.text(.02,.963,'Recall and precision in percent. Bold: best; italic: second-best in each predictor.',size=11,color='#536873')
        for i,predictor in enumerate(predictors):
            group=[r for r in rows if r['predictor']==predictor]
            ax=fig.add_axes([.02,.72-i*.235,.96,.208]);ax.axis('off')
            ax.text(0,1.04,NAMES[predictor],size=14,weight='bold',color='#266B8D',transform=ax.transAxes)
            values=[[NAMES[r['metric']]]+[r[c] for c in columns] for r in group]
            t=ax.table(cellText=values,colLabels=labels,cellLoc='center',colWidths=[.19]+[.09]*9,bbox=[0,0,1,1])
            t.auto_set_font_size(False);t.set_fontsize(10.2)
            for (j,k),cell in t.get_celld().items():
                cell.set_edgecolor('#D5DFE3');cell.set_linewidth(.35)
                cell.set_facecolor('#EAF0F3' if j==0 else '#F0F5F7' if group[j-1]['metric'] in ('nac','ric','gric') else 'white')
                if j==0:cell.get_text().set_weight('bold')
                if k==0:cell.get_text().set_ha('left')
            for k,col in enumerate(columns,1):
                ordered=sorted({float(r[col]) for r in group},reverse=col!='average_rank')
                for j,r in enumerate(group,1):
                    value=r[col]
                    if float(value)==ordered[0]:t[j,k].get_text().set_weight('bold')
                    elif len(ordered)>1 and float(value)==ordered[1]:
                        t[j,k].get_text().set_style('italic')
            ax.set_xlim(0,1);ax.set_ylim(0,1)
        fig.savefig(a.output/f'{stem}.png',dpi=170,bbox_inches='tight',facecolor='white')
        plt.close(fig)


if __name__=='__main__':main()
