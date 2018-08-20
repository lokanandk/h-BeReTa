function gTRE=h_BeReTa(model,BiomassRxn, ProductRxn, TRLevels, RegulatorGene, RegulatedGene,RegulatorySign, RegulatorGeneExp_Low,RegulatorGeneExp_High, RegulatedGeneExp_Low, RegulatedGeneExp_High, All_TRs)
%INPUTS 
%     Model                          Cobra model structure
%     BiomassRxn                     Biomass reaction
%     ProductRxn                     Product reaction
%     TRLevels                       TR hierarchies (Refer TR hierarchy level
%                                    information for E. coli provided as .mat file)
%     RegulatorGene                  String arrays of TRs 
%     RegulatedGene                  String arrays of the corresponding genes regulated
%     RegulatorySign                 +1 for Activator, -1 for Repressor, 0 for
%                                    Unknown/Dual regulator
%
%     RegulatorGeneExp_Low           Expression values of TRs for "Low or No"
%                                    product condition
%     RegulatorGeneExp_High          Expression values of TRs for "High"
%                                    product condition
%     RegulatedGeneExp_Low           Expression values of corresponging genes regulatedfor "Low or No"
%                                    product condition
%     RegulatedGeneExp_High          Expression values of corresponging genes regulatedfor "High"
%                                    product condition

% The sizes of all 4 quantities provided should be the same,
%i.e., equal to the total TR-gene interactions. Also, as explained in the
%main text and Figure 4, all loop forming interactions should be removed
%before providing them as inputs. For the case of E. coli, we have provided these TR-gene interactions for
%as .mat file, which can be used as inputs
nRAP=nRAPcalculation(model, BiomassRxn,ProductRxn);
[TRE,pvalues]=TREcalculation(RegulatorGene, RegulatedGene, RegulatorySign, RegulatorGeneExp_Low,RegulatorGeneExp_High, RegulatedGeneExp_Low,RegulatedGeneExp_High, nRAP);
gTRE=gTREcalculation(TRE,TRLevels,RegulatorGene, RegulatedGene, All_TRs,pvalues);
end


%%%Evaluation of nRAP scores
function nRAP = nRAPcalculation(model, BiomassRxn, ProductRxn)
%%Evaluation of nRAP scores
%%Decomposing reversible into two irreversible reactions
[modelIrrev] = convertToIrreversible(model);
modelIrrev.rxnNames=modelIrrev.rxns;
modelIrrev.metNames=modelIrrev.mets;
modelIrrev.description='Irrev_model';
modelIrrev = changeObjective(modelIrrev,BiomassRxn);
biomass_sol = optimizeCbModel(modelIrrev);
modelIrrev.lb(ismember(modelIrrev.rxns,BiomassRxn))=0.5*biomass_sol.f;
modelIrrev = changeObjective(modelIrrev, ProductRxn);
[minFlux,maxFlux]=fluxVariability(modelIrrev,50);
% %fastLooplessFVA avoid loops, but is much slower 
%[fvaRange] = fastLooplessFVA(modelIrrev,modelIrrev.c,0.5);
% %Reference: Saa, Pedro A., and Lars K. Nielsen. "Fast-SNP: a fast matrix pre-processing algorithm for efficient loopless flux optimization of metabolic models." Bioinformatics (2016): btw555.)
% %Here 0.5 indicates 50 % to optimality (alpha)
% minFlux=fvaRange(:,1);
% maxFlux=fvaRange(:,2); 
%%constrain-based formulation used by h-BeReTa
FluxAct2=zeros(length(modelIrrev.rxns),6);
for i=1:length(modelIrrev.rxns)    
        model1=modelIrrev;
        m=0;
        for k=0:0.2:1
        model1.lb(i)= minFlux(i)+ k*(maxFlux(i)-minFlux(i));
        model1.ub(i)=model1.lb(i);
        solution=optimizeCbModel(model1);
        m=m+1;
        if solution.f~=0
        FluxAct2(i,m)=solution.f;
        end
        end
end
%nRAP scoring for irreversible reactions
nRAP_ir=zeros(size(modelIrrev.rxns));
for j=1:length(modelIrrev.rxns)
    slope=gradient(FluxAct2(j,:),[0,0.2,0.4,0.6,0.8,1]);
    nRAP_ir(j,1)=slope(1);
end 
%%Averaging the nRAP scores of forward and backward reactions
nRAP=zeros(size(model.rxns));
m=1;
for l=1:(length(modelIrrev.rxns)-1)
    modelIrrev_new.rxns=strrep(modelIrrev.rxns,'_f','');
    modelIrrev_new.rxns=strrep(modelIrrev_new.rxns,'_b','');
    if strcmp(modelIrrev_new.rxns(l),modelIrrev_new.rxns(l+1))==1
        nRAP(m,1)=(nRAP_ir(l,1)+ nRAP_ir(l+1,1))/2;
    else
        if nRAP(m,1)==0
        nRAP(m,1)=nRAP_ir(l,1);
        end
        if m<length(model.rxns)
        m=m+1;
        end
        if l+1==length(modelIrrev.rxns)
            if nRAP(m,1)==0
                nRAP(m,1)=nRAP_ir(l+1);
            else
                nRAP(m,1)=(nRAP(m,1)+nRAP_ir(l+1))/2;
            end
        end
    end
    %%%%%%%%%% Note that the current code doesn't consist for nRAP evaluation for alternate flux modes
%%%%%%%%%% assisted by Fast-SL (Ref: )
%%%This part will be added soon
end
end


%%%TRE Evaluation
function [TRE,pvalues] = TREcalculation(RegulatorGene, RegulatedGene,RegulatorySign, RegulatorGeneExp_Low,RegulatorGeneExp_High, RegulatedGeneExp_Low,RegulatedGeneExp_High, nRAP)
nRS=zeros(length(RegulatedGene)); %normalized regulatory strengths
for g=1:length(RegulatedGene)
    nRS(g)=((RegulatedGeneExp_High - RegulatedGeneExp_Low)/(RegulatorGeneExp_High - RegulatorGeneExp_Low))*(RegulatorGeneExp_Low/RegulatedGeneExp_Low);
    if RegulatorySign(g)<0
        if nRS(g)>0
            nRS(g)=0;
        end
    end
    if RegulatorySign(g)>0
        if nRS(g)<0
            nRS(g)=0;
        end
    end
end
TRE=zeros(length(unique(RegulatorGene)),1);
OE=zeros(length(RegulatorGene),1);
for i=1:length(RegulatorGene)
    for j=1:length(nRAP)
        %%Note that in the current code can be used for E. coli iJO1366 model.
        %%Hence, GPRlevels and the corresponding GPR factors for this model
        %%are provided as .mat file, and hence need to be uploaded prior to
        %%using h-BeReTa
    if strcmp(RegulatedGene(i),GprPosition1(j))==1
        OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition2(j))==1
        OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition3(j))==1
        OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition4(j))==1
        OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition5(j))==1
        OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition6(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition7(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition8(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition9(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition10(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition11(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition12(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition13(j))==1
       OE(i)= OE(i) + ((nRAP(j)/GPRfactor(j))*nRS(i));
    end
    end
end
k=1;
for l=1:(length(RegulatorGene)-1)
if strcmp(RegulatorGene(l),RegulatorGene(l+1))==1
TRE(k,1)= TRE(k) + OE(l);
else
    TRE(k,1)= TRE(k) + OE(l);
    k=k+1;
end
end

%%%%%%Generating random nRAPs and corresponding TREs, this part consumes
%%%%%%significant amount of time (few hours) depending on the processor speed, code can be improved to optimize the speed 
RandTRE=zeros(length(unique(RegulatorGene)),1000);
RandOE=zeros(length(RegulatorGene),1000);
a=min(nRAP);
b=max(nRAP);
RandnRAP=a+(b-a)*rand(length(nRAP),1000);
for m=1:1000
for i=1:length(RegulatorGene)
    for j=1:length(nRAP)
    if strcmp(RegulatedGene(i),GprPosition1(j))==1
        RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition2(j))==1
        RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition3(j))==1
        RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition4(j))==1
        RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition5(j))==1
        RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition6(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition7(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition8(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition9(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition10(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition11(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition12(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    if strcmp(RegulatedGene(i),GprPosition13(j))==1
       RandOE(i,m)= RandOE(i,m) + ((RandnRAP(j,m)/GPRfactor(j))*nRS(i));
    end
    end
end
end
for n=1:1000
k=1;
for l=1:(length(RegulatorGene)-1)
if strcmp(RegulatorGene(l),RegulatorGene(l+1))==1
RandTRE(k,n)= RandTRE(k,n) + RandOE(l,n);
else
    RandTRE(k,n)= RandTRE(k,n) + RandOE(l,n);
    k=k+1;
end
end
end

%%% Calculating pvalues from actual and randomly generated TREs
pvalues=zeros(length(unique(RegulatorGene)),1);
for p=1:length(unique(RegulatorGene))
    if TRE(p,1)<0
    for q=1:1000
       pvalues(p,1)=length(find((RandTRE(p,:)<=(TRE(p,1)-0.15*TRE(p,1)))&(RandTRE(p,:)>=(TRE(p,1)+0.15*TRE(p,1)))))/1000;
    end
    elseif TRE(p,1)>0
    for q=1:1000
       pvalues(p,1)=length(find((RandTRE(p,:)>=(TRE(p,1)-0.15*TRE(p,1)))&(RandTRE(p,:)<=(TRE(p,1)+0.15*TRE(p,1)))))/1000;
    end
    end
end
end

%%%Global TRE calculation
function gTRE = gTREcalculation(TRE,TRLevels,RegulatorGene, RegulatedGene, All_TRs,pvalues)
TRs=unique(RegulatorGene);
TRE1=TRE;
gTRE=zeros(length(All_TRs));
for p=1:length(All_TRs)
    for q=1:length(TRs)
    if strcmp(All_TRs(p),TRs(q))==1
        if pvalues(q)<0.05
         gTRE(p)=TRE1(q);
        end
    end
    end
end
GbTempTRE=zeros(length(All_TRs),1);
for i=11:-1:1
for j=1:length(All_TRs)
for k=1:length(RegulatorGene)
if strcmp(TRLevels(j,i),RegulatorGene(k))==1
for m=1:length(All_TRs)
if strcmp(RegulatedGene(k),All_TRs(m))==1
for n=1:length(All_TRs)
if strcmp(RegulatorGene(k),All_TRs(n))==1
GbTempTRE(n,1)= nRS(k)*gTRE(m);
gTRE(n)=gTRE(n)+ GbTempTRE(n);   
end                    
end
end
end
end
end
end
end
end
%%%%%%%END of h-BeReTa