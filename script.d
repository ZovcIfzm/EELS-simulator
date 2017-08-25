import std.stdio, std.file, std.conv, arsd.dom, arsd.script, arsd.jsvar, main;
class Script
{
	Document file;
	PhaseSpace space;
	this(string fn)
	{
		this.file = new Document(readText(fn));
	}
	void run(){
		foreach(Element t; this.file.getElementsByTagName("space")){
			this.space = new PhaseSpace(to!double(t.getAttribute("width")), to!double(t.getAttribute("height")), 
				to!double(t.getAttribute("VzIntDist")), to!double(t.getAttribute("zIntDist")), 
				to!double(t.getAttribute("chirp")), to!double(t.getAttribute("b")), to!double(t.getAttribute("intensity")), to!double(t.getAttribute("intensityRatio")));
				this.byNode(t.childNodes()[0], this.space);
		}
	}
	void byNode(Element statement, PhaseSpace space){
		switch (statement.tagName) {
			case "model":
				space.modelPhaseSpace(to!double(statement.getAttribute("accuracy"))); 
				break;
			case "freeexpansion":
				space.freeExpansion(to!double(statement.getAttribute("time")));
				break;
			case "lens":
				space.opticalManipulation(to!double(statement.getAttribute("chirp")));
				break;
			case "print":
				space.printPhaseSpace();
				break;
			case "script":
				
				break;
			case "split":
				auto spaces = space.split(to!int(statement.getAttribute("portions")));
				foreach(PhaseSpace s; spaces){
					this.byNode(statement.childNodes()[0], s);
				}
				space = new PhaseSpace(spaces);
				break;
			default: break;
		}
		if(statement.nextSibling)
			this.byNode(statement.nextSibling, space);
	}
}
