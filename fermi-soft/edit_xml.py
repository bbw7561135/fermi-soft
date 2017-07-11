from xml.etree import ElementTree as ET

def run(name,P='1',I='0',S='0'):
    mod = ET.parse(str(name) + '_output_model.xml')
    model = mod.getroot()
    source = model.getchildren()

    #Find the source with the correct name
    i = 0
    while i < len(source):
        if model[i].attrib['name'] == str(name):
                params = source[i].getchildren()
        else:
            pass
        i += 1

    func = params[0]
    Pref = func[0]
    Index = func[1]
    Scale = func[2]

    Pref.attrib['free'] = P
    Index.attrib['free'] = I
    Scale.attrib['free'] = S

    mod.write(str(name) + '_null.xml')
