"""Reproduces the data figures of the paper (Figs. 2-7) from the CSV tables in data/tables.
Usage: python figures/make_figures.py   (run from the repository root; output: figures/out/*.pdf)"""
import os, math, csv
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, 'data', 'tables'); OUT = os.path.join(ROOT, 'figures', 'out') + os.sep
os.makedirs(OUT, exist_ok=True)
def load(name):
    with open(os.path.join(T, name)) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try: r[k] = float(v) if v not in ('', None) else None
            except ValueError: pass
    return rows
C=['#2a78d6','#eb6834','#1baf7a','#eda100']   # validated categorical slots 1-4
INK='#0b0b0b'; INK2='#52514e'; GRID='#e4e3df'
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.labelsize':8,
  'axes.titlesize':8,'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':0.6,
  'axes.edgecolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.labelcolor':INK,'text.color':INK,
  'xtick.major.width':0.6,'ytick.major.width':0.6,'lines.linewidth':1.4,'lines.markersize':4,
  'legend.frameon':False,'savefig.dpi':600,'pdf.fonttype':42})
def style(ax):
    ax.grid(True,which='major',color=GRID,lw=0.5); ax.set_axisbelow(True)
    for s in ('top','right'): ax.spines[s].set_visible(False)
def sup(x,pos=None):
    e=int(round(math.log10(x)));
    return r'$10^{%d}$'%e if abs(math.log10(x)-e)<1e-6 else ''
PL={'server':'Server','laptop':'Laptop'}
MK=['o','s','^','D']
def panel_label(ax,s): ax.text(0.5,-0.30,s,transform=ax.transAxes,ha='center',va='top',fontsize=8)

sw=load('bound_sweep.csv')
def rows(plat,b): return [r for r in sw if r['platform']==plat and r['bits']==b]

# ---------------- Fig: S0 scan (C++, 1280-bit)
s0=load('s0_scan_1280.csv')
fig,axs=plt.subplots(1,2,figsize=(7.16,2.05),sharey=False)
for ax,plat,wb in zip(axs,('server','laptop'),(64,32)):
    for j,r in enumerate((1.5,2,3)):
        d=[x for x in s0 if x['platform']==plat and float(x['r'])==r]
        bits=[x['first_product_bits'] for x in d]; y=[x['median_us'] for x in d]
        ax.plot(bits,y,color=C[j],marker=MK[j],ms=3,lw=1.2,label=f'$r={r:g}$')
    ax.axvline(wb,color=INK2,lw=0.8,ls=(0,(3,2)))
    ax.text(wb+2,ax.get_ylim()[1]*0.98 if False else None or max(y)*1.0,'',fontsize=7)
    ax.set_xlabel('Bit length of the first segment primorial'); ax.set_ylabel('Prefilter cost per candidate (µs)')
    style(ax); ax.legend(loc='upper left',ncol=3,handlelength=1.8)
    yl=ax.get_ylim(); ax.text(wb+1.5,yl[0]+0.04*(yl[1]-yl[0]),f'{wb}-bit word',color=INK2,fontsize=7)
panel_label(axs[0],'(a)'); panel_label(axs[1],'(b)')
fig.subplots_adjust(left=0.07,right=0.99,top=0.97,bottom=0.32,wspace=0.25)
fig.savefig(OUT+'fig_s0.pdf'); plt.close(fig)

# ---------------- Fig: B sweep: (a) prefilter cost 1280-bit; (b) TD/GEOM ratio four bit lengths (server)
fig,axs=plt.subplots(1,2,figsize=(7.16,2.2))
ax=axs[0]
for j,plat in enumerate(('server','laptop')):
    d=rows(plat,1280); B=[r['B'] for r in d]
    ax.plot(B,[r['segmented_prefilter_us'] for r in d],color=C[0],ls='-' if j==0 else (0,(4,2)),marker='o' if j==0 else 'o',mfc='white' if j else C[0],ms=3,label=f'Segmented, {PL[plat].lower()}')
    dt=[r for r in d if r['trial_division_us'] is not None]
    ax.plot([r['B'] for r in dt],[r['trial_division_us'] for r in dt],color=C[1],ls='-' if j==0 else (0,(4,2)),marker='s',mfc='white' if j else C[1],ms=3,label=f'Trial division, {PL[plat].lower()}')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Sieving bound $B$'); ax.set_ylabel('Prefilter cost per candidate (µs)')
style(ax); ax.legend(loc='upper left')
ax=axs[1]
for j,b in enumerate((768,1280,2560,3584)):
    for k,plat in enumerate(('server','laptop')):
        d=[r for r in rows(plat,b) if r['trial_division_us'] is not None]
        ax.plot([r['B'] for r in d],[r['trial_division_us']/r['segmented_prefilter_us'] for r in d],color=C[j],marker=MK[j],ms=3,
                ls='-' if k==0 else (0,(4,2)),mfc=C[j] if k==0 else 'white',label=f'{b} bits' if k==0 else None)
ax.axhline(1,color=INK2,lw=0.8)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Sieving bound $B$'); ax.set_ylabel('Cost ratio, trial division / segmented')
ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p: f'{v:g}')); ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([0.25,0.5,1,2,4,8])
style(ax); l1=ax.legend(loc='upper left',ncol=2)
from matplotlib.lines import Line2D
ax.add_artist(l1); ax.legend(handles=[Line2D([],[],color=INK2,ls='-',marker='o',ms=3),Line2D([],[],color=INK2,ls=(0,(4,2)),marker='o',mfc='white',ms=3)],labels=['Server','Laptop'],loc='lower right')
panel_label(axs[0],'(a)'); panel_label(axs[1],'(b)')
fig.subplots_adjust(left=0.07,right=0.99,top=0.97,bottom=0.28,wspace=0.27)
fig.savefig(OUT+'fig_sweep.pdf'); plt.close(fig)

# ---------------- Fig: expected cost relative to B=1000
fig,axs=plt.subplots(1,2,figsize=(7.16,2.1),sharey=True)
for ax,plat in zip(axs,('server','laptop')):
    for j,b in enumerate((768,1280,2560,3584)):
        d=rows(plat,b); cmr=d[0]['c_MR_us']
        B=np.array([r['B'] for r in d]); cost=np.array([r['segmented_prefilter_us']+r['pass_fraction_measured']*cmr for r in d])
        c0=cost[list(B).index(1000)]
        ax.plot(B,cost/c0,color=C[j],marker=MK[j],ms=3,label=f'{b} bits')
        i=int(np.argmin(cost)); ax.plot(B[i],cost[i]/c0,marker='*',ms=10,color=C[j],mec=INK,mew=0.5,zorder=5)
    ax.axhline(1,color=INK2,lw=0.8); ax.set_xscale('log'); ax.set_xlabel('Sieving bound $B$')
    style(ax)
axs[0].set_ylabel('Expected cost per candidate\nrelative to $B=1000$'); axs[0].legend(loc='upper center',ncol=2)
panel_label(axs[0],'(a)'); panel_label(axs[1],'(b)')
fig.subplots_adjust(left=0.08,right=0.99,top=0.97,bottom=0.31,wspace=0.08)
fig.savefig(OUT+'fig_expcost.pdf'); plt.close(fig)

# ---------------- Fig: B* scaling and formula accuracy
ob=load('best_bounds.csv'); op=load('best_bounds_prototype.csv')
fig,axs=plt.subplots(1,2,figsize=(7.16,2.2))
ax=axs[0]
sets=[('C++ library, server',[(r['bits'],r['Bstar_measured']) for r in ob if r['platform']=='server'],C[0],'o'),
      ('C++ library, laptop',[(r['bits'],r['Bstar_measured']) for r in ob if r['platform']=='laptop'],C[1],'s'),
      ('Prototype, laptop',[(r['bits'],r['Bstar_measured']) for r in op],C[2],'^')]
for name,pts,col,mk in sets:
    x=np.array([p[0] for p in pts],float); y=np.array([p[1] for p in pts],float)
    k,a=np.polyfit(np.log(x),np.log(y),1)
    xx=np.linspace(x.min(),x.max(),50); ax.plot(xx,np.exp(a)*xx**k,color=col,lw=1.0,ls=(0,(4,2)))
    ax.plot(x,y,mk,color=col,ms=4,label=f'{name} ($b^{{{k:.2f}}}$)')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Candidate bit length $b$'); ax.set_ylabel('Optimal sieving bound $B^*$')
ax.set_xticks([512,768,1024,1536,2048,3072,4096]); ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{int(v)}')); ax.xaxis.set_minor_formatter(NullFormatter())
style(ax); ax.legend(loc='upper left')
ax=axs[1]
xs=np.arange(8); lab=[f"{'S' if r['platform']=='server' else 'L'}{int(r['bits'])}" for r in ob]
ax.plot(xs-0.12,[r['closed_over_measured'] for r in ob],'o',color=C[0],ms=4.5,label='Closed form, $B\\log_2 B=c_{\\mathrm{MR}}/\\kappa$')
ax.plot(xs+0.12,[r['refined_over_measured'] for r in ob],'s',color=C[1],ms=4.5,label='Refined estimate')
for i,r in enumerate(ob): ax.plot([i-0.12,i+0.12],[r['closed_over_measured'],r['refined_over_measured']],color=GRID,lw=1,zorder=1)
ax.axhline(1,color=INK2,lw=0.8); ax.set_yscale('log'); ax.set_yticks([0.5,0.7,1,1.5,2,3]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}')); ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xticks(xs); ax.set_xticklabels(lab); ax.set_xlabel('Platform (S: server, L: laptop) and bit length'); ax.set_ylabel('Predicted / measured $B^*$')
style(ax); ax.legend(loc='upper right')
panel_label(axs[0],'(a)'); panel_label(axs[1],'(b)')
fig.subplots_adjust(left=0.08,right=0.99,top=0.97,bottom=0.28,wspace=0.25)
fig.savefig(OUT+'fig_bstar.pdf'); plt.close(fig)

# ---------------- Fig: safe primes
xs_=load('safe_prime_sweep.csv')
fig,axs=plt.subplots(1,2,figsize=(7.16,2.2))
ax=axs[0]
for j,b in enumerate((768,1280)):
    d=[r for r in xs_ if r['platform']=='server' and r['bits']==b and r['method']=='segmented_double_prefilter']
    B=np.array([r['B'] for r in d]); p=np.array([r['pass_fraction'] for r in d]); n=6000
    se=np.sqrt(p*(1-p)/n)
    ax.errorbar(B*(1.04 if j else 0.96),p,yerr=1.96*se,fmt=MK[j],color=C[j],ms=3.5,lw=0.8,capsize=1.5,label=f'Measured, {b} bits')
d=[r for r in xs_ if r['platform']=='server' and r['bits']==768 and r['method']=='segmented_double_prefilter']
B=np.array([r['B'] for r in d])
ax.plot(B,[r['W2_B'] for r in d],color=INK,lw=1.1,label=r'$W_2(B)\approx C_2W(B)^2$')
ax.plot(B,[r['W_B']**2 for r in d],color=INK2,lw=1.0,ls=(0,(4,2)),label=r'$W(B)^2$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Sieving bound $B$'); ax.set_ylabel('Fraction of pairs passing the prefilter')
ax.yaxis.set_minor_formatter(NullFormatter()); ax.set_yticks([0.004,0.006,0.01,0.02,0.03]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
style(ax); ax.legend(loc='upper right')
ax=axs[1]
for j,plat in enumerate(('server','laptop')):
    for k,m in enumerate(('segmented_double_prefilter','combined_sieve')):
        d=[r for r in xs_ if r['platform']==plat and r['bits']==1280 and r['method']==m]
        ax.plot([r['B'] for r in d],[r['expected_cost_us'] for r in d],color=C[k],marker=MK[k],ms=3,ls='-' if j==0 else (0,(4,2)),mfc=C[k] if j==0 else 'white',
                label=('Segmented double prefilter' if m=='segmented_double_prefilter' else 'Combined sieve')+f', {PL[plat].lower()}')
ax.set_xscale('log'); ax.set_xlabel('Sieving bound $B$'); ax.set_ylabel('Expected cost per candidate (µs)')
style(ax); ax.legend(loc='upper right')
panel_label(axs[0],'(a)'); panel_label(axs[1],'(b)')
fig.subplots_adjust(left=0.08,right=0.99,top=0.97,bottom=0.28,wspace=0.25)
fig.savefig(OUT+'fig_safe.pdf'); plt.close(fig)

# ---------------- Fig: confirmatory forest plot
cf=[r for r in load('openssl_confirmatory.csv') if r['primes_per_variant']]
# fixed-effect (pre-specified) and DerSimonian-Laird random-effects pooling of log(V0/V3)
yl=np.array([math.log(r['ratio_V0_over_V3']) for r in cf]); sl=np.array([r['se_log_ratio'] for r in cf]); w=1/sl**2
mfe=(w*yl).sum()/w.sum(); sfe=1/math.sqrt(w.sum()); Qc=(w*(yl-mfe)**2).sum()
tau2=max(0.0,(Qc-(len(yl)-1))/(w.sum()-(w**2).sum()/w.sum())); wr=1/(sl**2+tau2)
mre=(wr*yl).sum()/wr.sum(); sre=1/math.sqrt(wr.sum())
rows=[(r['ratio_V0_over_V3'],r['se_log_ratio'],False) for r in cf]+[(math.exp(mfe),sfe,True),(math.exp(mre),sre,True)]
lab=['Server, 1024 bits\n(OpenSSL 3.0.13)','Server, 1536 bits\n(OpenSSL 3.0.13)','Laptop, 1024 bits\n(OpenSSL 3.6.4)','Laptop, 1536 bits\n(OpenSSL 3.6.4)','Pooled, fixed effect\n(primary endpoint)','Pooled, random\neffects (sensitivity)']
fig,ax=plt.subplots(figsize=(3.5,2.6))
ys=np.arange(len(rows))[::-1]
for y,(m,s_,pooled),l in zip(ys,rows,lab):
    lo=math.exp(math.log(m)-1.96*s_); hi=math.exp(math.log(m)+1.96*s_)
    ax.plot([lo,hi],[y,y],color=INK2 if not pooled else INK,lw=1.2)
    ax.plot(m,y,'D' if pooled else 's',color=C[1] if pooled else C[0],ms=6 if pooled else 4.5,mec=INK if pooled else C[0],mew=0.5)
    ax.text(1.20,y,(f'{m:.3f} [{lo:.3f}, {hi:.3f}]' if pooled else f'{m:.3f} [{lo:.2f}, {hi:.2f}]'),va='center',fontsize=7)
ax.axvline(1,color=INK2,lw=0.8,ls=(0,(3,2)))
ax.set_yticks(ys); ax.set_yticklabels(lab,fontsize=7); ax.set_xlim(0.92,1.20); ax.set_xticks([0.95,1.0,1.05,1.10,1.15])
ax.set_xlabel('Ratio of mean generation times, stock / proposed')
style(ax); ax.grid(False,axis='y')
fig.subplots_adjust(left=0.30,right=0.73,top=0.98,bottom=0.15)
fig.savefig(OUT+'fig_confirm.pdf'); plt.close(fig)
print('ok')
