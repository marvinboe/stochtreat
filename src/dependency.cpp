/*
 *  dependency.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "dependency.h"

DependencyNode::DependencyNode(const DependencyNode& other){
    _affecting.clear();
    _idx = other.reaction();
    vector<DependencyNode*>::const_iterator start = other.begin();
    vector<DependencyNode*>::const_iterator stop = other.end();
    while(start !=stop){
        affects(*start);
        ++start;
    }
}

DependencyNode& DependencyNode::operator=(const DependencyNode& other){
    _affecting.clear();
    _idx =other.reaction();
    vector<DependencyNode*>::const_iterator start = other.begin();
    vector<DependencyNode*>::const_iterator stop = other.end();
    while(start !=stop){
        affects(*start);
        ++start;
    }
    return *this;
}

ostream& DependencyNode::display(ostream& os){
    os << "<" << _idx << "|"<< _affecting.size()<<"|";
    for(unsigned i=0; i< _affecting.size(); ++i){
        unsigned tmp =  _affecting[i]->reaction();
        os << tmp;
        if(i < (_affecting.size() -1))
            os << ",";
    }
    os << ">";
    return os;
}



//need to add a mutation reaction at a later stage; mutations occur between cell types. 
DependencyGraph::DependencyGraph(Model& pool, Data& data, AllReactions& all, RanGen& ran){
    //for every type of cell in every compartment there is one self-renewal and one differentiation reaction
    //note that compartment 0 is a special case, since it's size needs to remain constant
    unsigned pos;
    double sum=0.0;
    unsigned numtypes = 4;
    for(unsigned k=1; k < pool.numStoch(); ++k) { //every compartment

        //
        Treatment *treat= new Treatment(k,pool.getTreatRate());// type(k) -> type(k+1)+type(k+1)
        treat->setPropensity((treat->sufficientReactants(pool)?treat->reactantFactor(pool):0.0)); 
        sum += treat->propensity();
        pos=all.add(treat);
        DependencyNode *treatnode=new DependencyNode(pos);
        treatnode->affects(treatnode);
        pos=add(TREATMENT,treatnode);
        treat->setDG(TREATMENT,pos);

        for(unsigned type = 0; type < numtypes; ++type) {
            //			cout << "(1-eps) = " << (1.0-data.eps(type)) << "\t rate= " << pool.getRate(k) << endl; 
            SelfRenewal *self= new SelfRenewal(k,type,(1.0 - data.eps(type))*pool.getRate(k));// type(k) -> type(k)+type(k)
            self->setPropensity((self->sufficientReactants(pool)?self->reactantFactor(pool):0.0)); 
            sum += self->propensity();
            pos=all.add(self);
            DependencyNode *selfnode=new DependencyNode(pos);

            Differentation *diff= new Differentation(k,type,data.eps(type)*pool.getRate(k));// type(k) -> type(k+1)+type(k+1)
            diff->setPropensity((diff->sufficientReactants(pool)?diff->reactantFactor(pool):0.0)); 
            sum += diff->propensity();
            pos=all.add(diff);
            DependencyNode *diffnode=new DependencyNode(pos);


            //add dependencies between both nodes and themselves
            selfnode->affects(selfnode);
            selfnode->affects(diffnode);
            diffnode->affects(diffnode);
            diffnode->affects(selfnode);
            if (type==1) { // type is cancer cell or treated cell
                selfnode->affects(treatnode);
                diffnode->affects(treatnode);
            }
            if (type==1 || type==3){
                treatnode->affects(selfnode);
                treatnode->affects(diffnode);
            }


            //both nodes are also affected by certain reactions from the previous compartment.			
            //need to add them. Carefull with the reactions in compartment 0.
            if( (k-1) > 0){ //starting from k=2
                _diffnodes[((k-2)*4 + type)]->affects(selfnode);
                _diffnodes[((k-2)*4 + type)]->affects(diffnode);
                if (type==1){
                    _diffnodes[((k-2)*4 + type)]->affects(treatnode);
                }
            }
            //add both nodes to their corresponding vectors.
            pos=add(SELF_RENEWAL,selfnode);
            self->setDG(SELF_RENEWAL,pos);
            pos=add(DIFFERENTATION,diffnode);
            diff->setDG(DIFFERENTATION, pos);	
            //			cout << "Self renewal, K=" << k << ", reaction " << *self << endl;
            //			cout << "Differentiation, K=" << k << ", reaction " << *diff << endl;
        }
    }
    //add reactions for the stem cell compartment:
    for(int first_tp = 0; first_tp < 3; ++first_tp) { // excluding the bound type, i.e.tp = 3
        StemCellRenewal *scr=new StemCellRenewal(&ran,first_tp,pool.getRate(0));
        scr->setPropensity(scr->sufficientReactants(pool)?scr->reactantFactor(pool):0.0);
        sum += scr->propensity();
        pos=all.add(scr);
        DependencyNode *scrnode=new DependencyNode(pos);
        pos=add(MORAN,scrnode);
        scr->setDG(MORAN, pos);
        //		cout << "# stem cell renewal : " << *scr << "\t propensity " << <<endl;

        scrnode->affects(scrnode);
        //		cout << "#added loop dependency \n";

        if (_selfnodes.size() > 0 && _diffnodes.size() > 0){//TODO if first condition is true, other one automatically
            DependencyNode *nexttreat=get(TREATMENT,0); //TODO really needed for all types?
            scrnode->affects(nexttreat);
            for(int second_type = 0; second_type < 3; ++second_type) { //TODO why adding links to all cell types???
                DependencyNode *nextself, *nextdiff;
                nextself = get(SELF_RENEWAL,second_type); 
                nextdiff = get(DIFFERENTATION,second_type);
                scrnode->affects(nextself);
                //				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[nextself->reaction()]) << endl;
                scrnode->affects(nextdiff);
                //				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[nextdiff->reaction()]) << endl;
            }
            //			cout << "#added all dependencies between stem cell reactions and those in the next compartment\n";
        }

    }

    //add dependencies between all reactions in first (Moran) compartment
    for(unsigned first_tp = 0; first_tp < 3; ++first_tp) { 
        DependencyNode *scrnode =  get(MORAN,first_tp);
        for(unsigned second_tp = 0; second_tp < 3; ++second_tp) {
            if(first_tp != second_tp){
                DependencyNode *oscrnode =  get(MORAN,second_tp);
                scrnode->affects(oscrnode);
                //				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[oscrnode->reaction()]) << endl;

            }
        }
    }
    all.setPropSum(sum);
    // cout << "#finished creation DependencyGraph, number of reactions = "<< all.size() << endl;
}


ostream& DependencyGraph::display(ostream& os){
    os << "DependencyGraph [" << endl << "Moran reactions";
    vector<DependencyNode*>::iterator startmoran = _morannodes.begin();
    vector<DependencyNode*>::iterator stopmoran = _morannodes.end();
    while(startmoran !=stopmoran){
        os << endl << **startmoran;
        ++startmoran;
    }
    os << endl<< "Selfrenewal";
    vector<DependencyNode*>::iterator startself = _selfnodes.begin();
    vector<DependencyNode*>::iterator stopself = _selfnodes.end();
    while(startself !=stopself){
        os << endl << **startself;
        ++startself;
    }
    os << endl << "Differentiation";
    vector<DependencyNode*>::iterator startdiff = _diffnodes.begin();
    vector<DependencyNode*>::iterator stopdiff = _diffnodes.end();
    while(startdiff !=stopdiff){
        os << endl << **startdiff;
        ++startdiff;
    }
    os << "]";
    return os;
}




