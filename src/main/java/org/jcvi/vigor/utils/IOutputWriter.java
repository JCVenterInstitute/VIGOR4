package org.jcvi.vigor.utils;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;

import java.io.IOException;
import java.util.List;

public interface IOutputWriter {

    void writeModels(Outfiles outfiles, List<Model> models) throws VigorException, IOException;
    String getExtension();
    WriterBundle getWriter(Outfiles outfiles,
                           OutputContext context,
                           OutputContext.Key ... toClose ) throws IOException, VigorException;
}
