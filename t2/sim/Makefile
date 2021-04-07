all: ngspice1.log ngspice2.log ngspice3.log ngspice4.log ngspice5.log
	

	$(eval NGSTABS1=$(shell grep _TAB ngspice1.log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS1), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' ngspice1.log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	
	
	$(eval NGSTABS2=$(shell grep _TAB ngspice2.log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS2), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' ngspice2.log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)


	$(eval NGSFIGS3=$(shell grep _FIG ngspice3.log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS3=$(addsuffix .ps, $(NGSFIGS3)))
	$(foreach i, $(NGSFIGPS3), ps2pdf $i;)
	

	$(eval NGSFIGS4=$(shell grep _FIG ngspice4.log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS4=$(addsuffix .ps, $(NGSFIGS4)))
	$(foreach i, $(NGSFIGPS4), ps2pdf $i;)
	
	
	$(eval NGSFIGS5=$(shell grep _FIG ngspice5.log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS5=$(addsuffix .ps, $(NGSFIGS5)))
	$(foreach i, $(NGSFIGPS5), ps2pdf $i;)


ngspice1.log: ../sim/sim1.net
	ngspice -b $< -o $@

ngspice2.log: ../sim/sim2.net                         
	ngspice -b $< -o $@

ngspice3.log: ../sim/sim3.net                          
	ngspice -b $< -o $@

ngspice4.log: ../sim/sim4.net                        
	ngspice -b $< -o $@
	
ngspice5.log: ../sim/sim5.net                        
	ngspice -b $< -o $@


clean:
	@rm -f *.log *.tex *.pdf *.ps *.txt *~

.PHONY: all clean
