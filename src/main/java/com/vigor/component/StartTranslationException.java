package com.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import java.util.List;

/**
 * Created by snettem on 5/8/2017.
 */

@Component
@Scope("prototype")
@Data
public class StartTranslationException {

   private boolean isStartTranslationException=false;
   private List<String> alternateStartCodons;

}
