import json
from sys import exit

#bporder = ('19','20','28','11','08','23','16','24','03','06','07','10','18',
#           '02','05','04','12','15','21','01','90','25','29','13','27')

#====== Residue ================================================================
class Residue:
    pass

#====== Mplets =================================================================
class Mplet:
    pass

#====== Stem ===================================================================
class Stem:
    pass
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def revbps(self):
        return [(nt2,nt1) for nt1,nt2 in self.bps[::-1]]

#====== Helices ================================================================
class Helix:
    pass
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def revbps(self):
        return [(nt2,nt1) for nt1,nt2 in self.bps[::-1]]

#====== Hairpins ===============================================================
class Hairpin:
    name = 'hairpin'

#====== Bulges =================================================================
class Bulge:
    name = 'bulge'

#====== Iloops =================================================================
class Iloop:
    name = 'iloop'

#====== Junction ===============================================================
class Junct:
    name = 'junct'

#====== NewJunction ============================================================
class Newj:
    name = 'newj'

#====== Chain ==================================================================
class Chain:
    pass

#===============================================================================
class Dssr:
    pass

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def __init__(self, fn=None):
        self.nts = {}
        self.bps = []
        self.bpinfo = {}
        self.mplets = []
        self.stems = {}
        self.stacks = []
        self.helices = {}
        self.hps = []
        self.bulges = []
        self.iloops = []
        self.juncts = []
        self.newjs = []
        self.chains = []
        if fn:
            self.read(fn)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def read(self, fn):
        # ----------------------------------------------------------------------
        def foundfirst(loop, nts):
            # useful when searching for non-canonical bps in junction loops
            # Purpose: find the first occurence of a loop nt that is included
            #          in nts.
            found = False
            for i,nt in enumerate(loop):
                if nt in nts:
                    found = True
                    idx = i
                    break
            if found:
                result = idx
            else:
                result = -1
            return result
        # ----------------------------------------------------------------------
        #
        f = open(fn)
        data = json.load(f)
        f.close()
        # if some keys in data is missing, add it with an empty value
        for k in ('nts','pairs','multiplets','helices','stems','coaxStacks',
                  'stacks','isoCanonPairs','hairpins','bulges','iloops',
                  'junctions'):
            if k not in data:
                data[k] = []
        # construct residues
        attrs = ('alpha','beta','gamma','delta','epsilon','zeta','chi')
        for nt in data['nts']:
            res = Residue()
            res.chain = nt['chain_name']
            res.resi = nt['nt_resnum']
            res.name = nt['nt_name']
            res.form = nt['form']
            res.splay = (nt['splay_angle'], nt['splay_distance'],
                         nt['splay_ratio'])
            res.phase = nt['phase_angle']
            res.pucker = nt['puckering']
            res.sugar = nt['sugar_class']
            # to make it work for dssr v1.X
            try:
                res.type = nt['nt_type']
            except KeyError:
                res.type = None
            # add attr in the attrs list
            for attr in attrs:
                setattr(res, attr, nt[attr])
            try:
                res.glyco = nt['glyco_bond']
            except KeyError:
                res.glyco = nt['summary'].split(',')[0]
            self.nts[nt['nt_id']] = res
        # construct base pairs
        for pair in data['pairs']:
            nt1,nt2 = pair['nt1'],pair['nt2']
            self.bps.append((nt1, nt2))
            info = [pair['bp'], pair['name'], pair['Saenger'][:2], pair['LW']]
            self.bpinfo[(nt1,nt2)] = info
            self.bpinfo[(nt2,nt1)] = info
        # construct multiplets
        for mp in data['multiplets']:
            mplet = Mplet()
            self.mplets.append(mplet)
            mplet.seq = tuple(mp['nts_long'].split(','))
            mplet.s = mp['nts_short']
            # to make it work for dssr v1.X
            try:
                mplet.planar = mp['planarity']
            except KeyError:
                mplet.planar = None
        # construct helices
        for hx in data['helices']:
            helix = Helix()
            helix.sidx = []
            helix.iidx = []
            idx = hx['index']
            self.helices[idx] = helix
            helix.s1 = hx['strand1']
            helix.s2 = hx['strand2']
            helix.type = hx['bp_type']
            helix.hform = hx['helix_form']
            helix.bps = [(pair['nt1'],pair['nt2']) for pair in hx['pairs']]
        # construct stems
        for sm in data['stems']:
            stem = Stem()
            idx = sm['index']
            self.stems[idx] = stem
            stem.s1 = sm['strand1']
            stem.s2 = sm['strand2']
            stem.type = sm['bp_type']
            stem.hform = sm['helix_form']
            stem.bps = [(pair['nt1'],pair['nt2']) for pair in sm['pairs']]
            stem.hidx = sm['helix_index']
            # add stem index to the corresponding helix object
            self.helices[stem.hidx].sidx.append(idx)
        # update coaxial-stacking info in helice[].sidx
        for ca in data['coaxStacks']:
            # ----------> for sanity check, and can be removed later -----------
            if len(self.helices[ca['helix_index']].sidx) != \
               len(ca['stem_indices']):
                print('Design ERROR! Check coaxStacks parsing code.')
                exit(1)
            # <-----------------------------------------------------------------
            self.helices[ca['helix_index']].sidx = ca['stem_indices']
        # construct stacks
        for sk in data['stacks']:
            self.stacks.append(tuple(sk['nts_long'].split(',')))
        # construct isobps and merge it into stem with negative idx
        for ibp in data['isoCanonPairs']:
            stem = Stem()
            idx = ibp['index']
            self.stems[idx] = stem
            stem.bps = [(ibp['nt1'], ibp['nt2'])]
            stem.hidx = ibp['helix_index']
            # when isobp does not belong to a helix, hidx will be out of range
            # In this case, assign None to stem.hidx
            if stem.hidx not in self.helices.keys():
                stem.hidx = None
                print('An isobp does not belong to a helix. Please check it!!!')
                exit(1)
            # add negative stem index to the corresponding helix object
            self.helices[stem.hidx].iidx.append(idx)
        # construct hairpins
        for hp in data['hairpins']:
            hpin = Hairpin()
            self.hps.append(hpin)
            hpin.sidx = hp['stem_indices'][0]
            hpin.s = hp['bridges'][0]['nts_short']
            hpin.seq = tuple(hp['bridges'][0]['nts_long'].split(','))
            hpin.nts = hp['nts_long'].split(',')
            hpin.closing = (hpin.nts[0], hpin.nts[-1])
            hpin.stem = self.stems[hpin.sidx].bps
            # ----------> for sanity check, and can be removed later -----------
            if len(hp['bridges'])!=1 or len(hp['stem_indices'])!=1:
                print('Design ERROR! Check hairpin parsing code.')
                exit(1)
            # <-----------------------------------------------------------------
        # construct bulges
        for bg in data['bulges']:
            bulge = Bulge()
            self.bulges.append(bulge)
            left = bg['bridging_nts'][0]
            loop = [x for x in bg['bridges'] if x['num_nts']!=0][0]
            bulge.sidx = bg['stem_indices']
            bulge.s = loop['nts_short']
            bulge.seq = tuple(loop['nts_long'].split(','))
            nts = tuple(bg['nts_long'].split(','))
            bulge.stems = []
            if left:
                # if bulge is on the left
                bulge.nts = nts
                bulge.closing = [(nts[0],nts[-1]), (nts[-3],nts[-2])]
                for si in bulge.sidx:
                    bulge.stems.append(self.stems[si].bps)
            else:
                # if bulge is on the right, then convert it to the left
                bulge.nts = tuple(list(nts[2:]) + list(nts[:2]))
                bulge.closing = [(nts[2],nts[1]), (nts[-1],nts[0])]
                for si in bulge.sidx[::-1]:
                    bulge.stems.append(self.stems[si].revbps())
        # construct internal loops
        for il in data['iloops']:
            iloop = Iloop()
            self.iloops.append(iloop)
            n1 = il['bridging_nts'][0]
            iloop.sidx = il['stem_indices']
            iloop.s = [x['nts_short'] for x in il['bridges']]
            iloop.seq = [tuple(x['nts_long'].split(',')) for x in il['bridges']]
            nts = tuple(il['nts_long'].split(','))
            iloop.nts = nts
            iloop.closing = [(nts[0],nts[-1]), tuple(nts[n1+1:n1+3])]
            iloop.stems = []
            for si in iloop.sidx:
                iloop.stems.append(self.stems[si].bps)
        # construct junctions
        for jt in data['junctions']:
            junct = Junct()
            self.juncts.append(junct)
            nway = int(jt['type'].split('-')[0])
            junct.sidx = jt['stem_indices']
            # reverse bp order of 1st stem
            junct.stems = [self.stems[junct.sidx[0]].revbps()]
            # append the rest of stems
            junct.stems += [self.stems[i].bps for i in junct.sidx[1:]]
            # retrieve linkers
            junct.loops = [tuple(b['nts_long'].split(','))
                           for b in jt['bridges']]
        # construct my new junctions based on those from dssr
        for junct in self.juncts:
            newj = Newj()
            self.newjs.append(newj)
            newj.sidx = junct.sidx
            newj.stems = junct.stems  # stems remain unchanged
            newj.exts = []
            newj.bulges = []
            loops = list(junct.loops)  # loops are modified from old loops
            nstem = len(newj.stems)
            for i in range(nstem):
                floop = loops[i]  # forward loop
                bloop = loops[i-1+nstem if i-1<0 else i-1]  # backward loop
                si = newj.sidx[i]  # fetch stem#
                hi = self.stems[si].hidx  # fetch helix# that contains this stem
                H = self.helices[hi].bps
                if junct.stems[i][0] not in H:
                    H = self.helices[hi].revbps()
                Hnts = [nt for bp in H for nt in bp]
                # gradually search for bp in floop + bloop
                ext = []
                bulge = False
                while True:
                    loc1 = foundfirst(floop, Hnts)
                    if loc1 == -1:
                        break
                    loc2 = foundfirst(bloop[::-1], Hnts)
                    if loc2 == -1:
                        break
                    if loc1 or loc2:  # loc1 or loc2 nonzero: bulges exist
                        bulge = True
                    bp = (bloop[::-1][loc2], floop[loc1])
                    if bp not in H:
                        break
                    ext.append(bp)
                    floop = floop[loc1+1:]
                    bloop = bloop[::-1][loc2+1:][::-1]
                loops[i] = floop
                loops[i-1+nstem if i-1<0 else i-1] = bloop
                newj.exts.append(ext)
                newj.bulges.append(bulge)
            newj.loops = loops
        # construct chains
        for k in data['chains']:
            chain = Chain()
            self.chains.append(chain)
            ch = data['chains'][k]
            chain.name = k.split('_')[-1]
            chain.s = ch['bseq']
            chain.dbn = ch['sstr']
            chain.form = ch['form']

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def pairing(self, nt):
        return [bp for bp in self.bps if nt in bp], \
               [mp.seq for mp in self.mplets if nt in mp.seq]

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def searchstem(self, sidx):
        # Given a stem#, return motifs containing this stem
        result = []
        for hp in self.hps:
            if sidx == hp.sidx:
                result.append(hp)
        for bg in self.bulges:
            if sidx in bg.sidx:
                result.append(bg)
        for il in self.iloops:
            if sidx in il.sidx:
                result.append(il)
        for jt in self.juncts:
            if sidx in jt.sidx:
                result.append(jt)
        return result

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
